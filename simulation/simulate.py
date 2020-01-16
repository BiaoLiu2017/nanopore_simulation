from simulation.SimSignalGenerator import SimSignalGenerator
import numpy as np
import h5py
import random
from datetime import date
import time
import os


def run(parser, args):
    my_generator = SimSignalGenerator(config_file=args.config)
    my_generator.load_reference(ref_file=args.ref)
    my_generator.load_model(args.model)

    start_time = int(time.time())
    directory = args.output
    if not os.path.exists(directory):
        os.makedirs(directory)
    (N, n) = divmod(args.reads,args.dirreads)
    for k in range(N + 1):
        if k == N:
            if n == 0:
                break
            number_reads = n
        else:
            number_reads = args.dirreads
        new_directory = os.path.join(directory, str(k))
        if not os.path.exists(new_directory):
            os.makedirs(new_directory)
        j = 0
        (M, m) = divmod(number_reads,args.workerreads)

        for j in range(M+1):
            if j == M:
                sub_number_reads = m
            else:
                sub_number_reads = args.workerreads

            reads = my_generator.generate(sub_number_reads, sampling_rate=my_generator.signal_config.sampling_rate, \
                cut_off_freq=args.cut_off, bandwidth_freq=args.band, error_rate=args.error_rate, is_reverse=args.reverse, \
                    is_signal_repeat=args.signal_repeat, is_uniform_noise=args.uniform_noise, is_event_repeat=args.event_repeat, \
                        is_low_pass_filter=args.low_pass_filter,is_gaussian_noise=args.gaussian_noise, is_correction=args.correction)
            random.shuffle(reads)

            for i, read in enumerate(reads):
                values_offset = read[2]
                values_range = read[3]
                sampling_rate = my_generator.signal_config.sampling_rate
                signal = read[1]
                read_number = k*args.dirreads + j*args.workerreads + i + 1
                hostname = "Simulator"
                date_string = date.today().strftime("%Y%m%d")
                device = "Nanopore_SimulatION"
                filename = hostname + "_" + date_string + "_" + device + "_" + "_read" + str(read_number) + "_strand.fast5"

                f = h5py.File(os.path.join(new_directory, filename), 'w')

                grp = f.create_group("UniqueGlobalKey/channel_id")
                grp.attrs.create("offset", data=values_offset, dtype="float64")
                grp.attrs.create("range", data=values_range, dtype="float64")
                grp.attrs.create("sampling_rate", data=my_generator.signal_config.sampling_rate, dtype="float64")

                grp = f.create_group("UniqueGlobalKey/context_args")
                grp.attrs.create("filename", data=np.string_(filename))
                grp.attrs.create("sample_frequency", data=np.string_(str(my_generator.signal_config.sampling_rate)))

                grp = f.create_group("UniqueGlobalKey/tracking_id")
                grp.attrs.create("device_id", data=np.string_(device))
                grp.attrs.create("exp_start_time", data=np.string_(str(start_time)))
                grp.attrs.create("hostname", data=np.string_(hostname))

                grp = f.create_group("Raw/Reads/Read_" + str(read_number))
                grp.attrs.create("sequence", data=np.string_(read[0]))
                grp.attrs.create("duration", data=len(signal), dtype="int32")
                grp.attrs.create("median_before", data=250, dtype="float64")
                grp.attrs.create("read_id",
                                 data=np.string_("16acf7fb-696b-4b96-b95b-0a43f" + format(read_number, '07d')))
                grp.attrs.create("read_number", data=read_number, dtype="int32")
                grp.attrs.create("start_mux", data=1, dtype="int32")
                grp.attrs.create("start_time", data=int(time.time()), dtype="int64")

                grp.create_dataset("Signal", data=signal.astype("int16"), dtype="int16", compression="gzip",
                                   compression_opts=1, maxshape=(None,))

                f.close()
            print("Simulated " + str(read_number) + " of " + str(args.reads) + " reads so far...")