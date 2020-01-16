import pickle
import numpy as np
def read_pkl(pkl_file):
    with open(pkl_file, "rb") as fp:   # Unpickling
        b = pickle.load(fp)
        return b

def write_pkl(contents, pkl_file):
    with open(pkl_file, "wb") as fp:   #Pickling
        pickle.dump(contents, fp)

def main():
    pklfile = "/home/liubiao/Chiron/Simulator/nanopore_simulation/examples/DNA/gDNAWA01.sic"
    contents = read_pkl(pklfile)
    '''
    6:Event length distribution
    7:offsets
    8:ranges
    9:Sampling rate
    '''
    contents = np.array(contents)[[6,7,8,9]].tolist()
    config_file = "/home/liubiao/Chiron/Simulator/nanopore_simulation/examples/DNA/config.txt"
    write_pkl(contents, config_file)
    print(read_pkl(config_file))

if __name__ == '__main__':
    main()