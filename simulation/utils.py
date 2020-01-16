from Bio import SeqIO
from Bio import Seq
import numpy as np
from random import *
import scipy.signal

def reverse(seq):
    seq = Seq.reverse_complement(seq)
    return seq

def mutation(new_record, amount_error, base_list):
    mis_error, ins_error, del_error = 0, 0, 0
    error_sites = np.random.choice(len(new_record), amount_error, replace=False)
    error_sites.sort()
    for error_site in error_sites:
        # Random error type 1 = Mismatch, 2 = Insertion, 3 = Deletion
        error_type = randint(1, 3)
        if error_type == 1:
            mis_error += 1
            choice = sample(base_list, 2)
            if new_record[error_site] == choice[0]:
                new_record[error_site] = choice[1]
            else:
                new_record[error_site] = choice[0]

        elif error_type == 2:
            ins_error += 1
            choice = sample(base_list, 2)
            new_record.insert(error_site,choice[0])
            error_sites[error_sites > error_site] = error_sites[error_sites > error_site] + 1
        else:
            del_error += 1
            del(new_record[error_site])
            error_sites[error_sites > error_site] = error_sites[error_sites > error_site] - 1
    return new_record, mis_error, ins_error, del_error

def signal_repeat(kmer_means):
    event_move = np.ceil(np.random.exponential(scale=2, size=len(kmer_means))).astype('int')
    event_idx = np.repeat(np.arange(len(kmer_means)), event_move)
    event_mean = kmer_means[event_idx]
    return event_mean, event_idx

def event_repeat(event_mean, event_samples, event_idx):
    event_move = event_samples
    event_idx = np.repeat(event_idx, event_move)
    event_repeat_idx = np.repeat(np.arange(len(event_mean)), event_move)
    event_mean = event_mean[event_repeat_idx]
    return event_mean, event_idx

def uniform_noise(kmer_stdvs, event_idx):
    uniform_noise = np.random.uniform(-1*kmer_stdvs[event_idx], kmer_stdvs[event_idx])
    return uniform_noise

def gaussian_noise(event_mean):
    gaussian_noise = np.random.normal(0, 1, len(event_mean))
    return gaussian_noise

def low_pass_filter(sampling_rate, cut_off_freq, bandwidth_freq):
    # Read input parameter
    fS = sampling_rate  # Sampling rate.
    fL = cut_off_freq   # Cutoff frequency.
    fb = bandwidth_freq # Bandwidth frequency

    # Generate frequency bin
    b = fb / fS
    N = int(np.ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = np.arange(N)

    # Compute sinc filter.
    h = np.sinc(2 * fL / fS * (n - (N - 1) / 2.))

    # Compute Blackman window.
    w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + \
        0.08 * np.cos(4 * np.pi * n / (N - 1))

    # Compute h and h_start
    h = h * w
    h /= np.sum(h)
    impulse = np.repeat(0., len(h))
    impulse[0] = 1.
    h_response = scipy.signal.lfilter(h, 1, impulse)
    h_start = np.argmax(h_response)

    # return
    return h,h_start,N