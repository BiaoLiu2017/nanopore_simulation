from .SIConfigFile import SIConfigFile
from simulation import utils
from random import *
from Bio import SeqIO
from Bio import Seq
import numpy as np
from pylab import *
import scipy.signal

class SimSignalGenerator(object):

    def __init__(self,
                 config_file='',
                 debug=False):
        self.new_record = ""
        self.debug = debug
        self.ref_data = ""
        self.ref_data_keys = []
        self.ref_length = 0
        self.base = []
        self.k = 0
        self.case = 0
        self.min_length = 0
        self.max_length = 0
        self.signal_config = None
        self.read_length_distribution = []
        self.event_length_distribution = []
        self.offsets = []
        self.ranges = []
        self.bases_per_second = 0
        self.pores_number = 0
        self.max_active_pores = 0
        self.read_until = None
        self.wear_out = []
        self.file_path = None
        self.empty_regions = False
        self.model_data = None
        self.model_dict = {}
        self.load_config(config_file)

    def load_model(self, model_file):
        print("Loading model file")
        self.model_data = np.genfromtxt(model_file, delimiter="\t", dtype=None, comments='#', names=True)
        # k = kmer length
        self.k = len(self.model_data[0][0])
        self.model_dict = dict([(x[0].decode('utf-8'), (x[1], x[2])) for x in self.model_data])

    def load_reference(self, ref_file, low_mem=False):
        print("Loading reference genome. May take some time depending on genome size.")
        if low_mem:
            self.ref_data = SeqIO.index(ref_file, "fasta")
        else:
            self.ref_data = SeqIO.to_dict(SeqIO.parse(ref_file,"fasta"))
        self.case = 0
        # Base dictionary for choosing error bases
        self.base = ['A', 'C', 'G', 'T']
        if len(self.ref_data) > 1:
            print("More than one sequence, but only simulate the first sequence!")
        elif len(self.ref_data) == 0:
            print("There is no sequence in the reference file!")
        for key in self.ref_data:
            self.ref_data_keys.append(key)
            # if self.ref_data[key].seq.count("M") > 0:
            #     if self.debug: print("Bases M detected in sequence")

            #     # Set case to 1
            #     self.case = 1
            #     # New Base Dictionary in case 1 with extra bases M, N for choosing error bases
            #     self.base = ['A', 'C', 'G', 'T', 'M']
            if self.ref_data[key].seq.count("N") > 0:
                self.empty_regions = True

    def load_config(self, config_file):
        print("Loading config file")
        handle = open(config_file,mode='rb')
        self.signal_config = SIConfigFile()
        self.signal_config.load_file(handle)
        # lengths = []
        # probabilities = []
        # for rld in self.signal_config.read_length_distribution:
        #     lengths.append(rld[0])
        #     probabilities.append(rld[1])
        # self.read_length_distribution.append(lengths)
        # self.read_length_distribution.append(probabilities)
        lengths = []
        probabilities = []
        for eld in self.signal_config.event_length_distribution:
            lengths.append(eld[0])
            probabilities.append(eld[1])
        self.event_length_distribution.append(lengths)
        self.event_length_distribution.append(probabilities)
        lengths = []
        probabilities = []
        for offset in self.signal_config.offsets:
            lengths.append(offset[0])
            probabilities.append(offset[1])
        self.offsets.append(lengths)
        self.offsets.append(probabilities)
        lengths = []
        probabilities = []
        for range in self.signal_config.ranges:
            lengths.append(range[0])
            probabilities.append(range[1])
        self.ranges.append(lengths)
        self.ranges.append(probabilities)
        # self.bases_per_second = self.signal_config.bases_per_second
        # self.pores_number = self.signal_config.pores_number
        # self.max_active_pores = self.signal_config.max_active_pores
        # self.read_until = self.signal_config.read_until
        # self.wear_out = self.signal_config.wear_out
        # self.file_path = self.signal_config.file_path


    def generate(self, snip_count=1, precise=False, debug=False, sampling_rate=4000.0, cut_off_freq=1750.0, 
    bandwidth_freq=40.0, error_rate=0.0, is_reverse=False, is_signal_repeat=True, is_uniform_noise=True, is_event_repeat=True, 
        is_low_pass_filter=True, is_gaussian_noise=True, is_correction=True):
        # Generate snippets
        first_key = next(iter(self.ref_data))
        snippets = np.array([self.ref_data[first_key].seq])
        reads = []
        for i in range(snip_count):
            seq = snippets[0]
            if is_reverse:
                seq = ''.join(list(seq))
                seq = utils.reverse(seq)
            new_record = list(seq)

            # Generate Errors
            amount_error = int(error_rate * len(new_record))
            if amount_error > 0:
                new_record, mis_error, ins_error, del_error = utils.mutation(new_record=new_record, amount_error=amount_error, base_list=self.base)
                print("Amount mismatch: %s, insertion: %s, deletion: %s" % (mis_error, ins_error, del_error))
            
            # New Sequence
            sequence = ""
            # if in the reference not all bases are known these need to be simulated as well
            if self.empty_regions:
                new_record = np.array(new_record)
                new_record[new_record == 'N'] = np.random.choice(self.base, 1)
            sequence = "".join(new_record)

            # Divide sequence into kmers
            kmers = [sequence[i:i + self.k] for i in range(0, len(sequence) - self.k + 1)]
            try:
                kmer_means, kmer_stdvs = zip(*[self.model_dict[kmer] for kmer in kmers])
            except ValueError:
                print("Reference too short to simulate. Continueing!")
                continue
            kmer_means = np.array(kmer_means)
            kmer_stdvs = np.array(kmer_stdvs)

            #signal repeat
            if is_signal_repeat:
                event_mean, event_idx = utils.signal_repeat(kmer_means)
            else:
                event_mean = kmer_means
                event_idx = np.arange(len(kmer_means))
            
            #uniform noise
            if is_uniform_noise:
                event_std = utils.uniform_noise(kmer_stdvs, event_idx)
                event_mean = event_mean + event_std
            else:
                event_std = utils.uniform_noise(kmer_stdvs, event_idx)

            #event repeat
            if is_event_repeat:
                event_total = len(event_mean)
                event_samples = np.random.choice(self.event_length_distribution[0], event_total, p=self.event_length_distribution[1]).astype(int)
                event_mean, event_idx = utils.event_repeat(event_mean, event_samples, event_idx)
            else:
                event_samples = np.array((1, )*len(event_mean))
            
            #low pass filter
            if is_low_pass_filter:
                h,h_start,N = utils.low_pass_filter(sampling_rate, cut_off_freq, bandwidth_freq)
                event_mean = np.convolve(event_mean.tolist(),h)[h_start+1:-(N-h_start-1)+1]

            #gaussian noise
            if is_gaussian_noise:
                event_list = np.stack([np.abs(event_std), event_samples], axis=1)
                Noise = []
                [Noise.extend(np.random.normal(((float(event[0])) * -1/2), (2*float(event[0])/2), int(event[1]))) for event in event_list]
                event_mean = event_mean + Noise
            else:
                gaussian_noise = utils.gaussian_noise(event_mean)
                event_mean = event_mean + gaussian_noise
            
            #correction
            if is_correction:
                signal_offset = float(np.random.choice(self.offsets[0],1,p=self.offsets[1]))
                signal_range = np.random.choice(self.ranges[0],1,p=self.ranges[1])
                event_mean = np.ceil(event_mean * (8192 / signal_range) + signal_offset)
            
            signal = event_mean

            reads.append([sequence, np.array(signal).astype(int), float(signal_offset), float(signal_range)])

        return reads
