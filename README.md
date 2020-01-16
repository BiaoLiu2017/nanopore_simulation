## Nanopore SimulatION

| Inputs | Outputs |
|---|---|
| Configuration-File (derived form a real conducted experiment)| 
| Reference genome FASTA-File (may be from another species)| Basecallable Fast5-Files with raw values |
| Model-File (provided by ONT)| |

### Description

Nanopore SimulatION is a tool for simulating an Oxford Nanopore Technologies MinION device for bioinformatic development. 

### Installation

Install the Nanopore SimulatION package:
```bash
conda create -n Simulator python=3
conda activate Simulator
```

```bash
git clone ssh://git@192.168.224.185:7999/~liubiao2/simulator.git
cd simulator
pip install -e ./
```

##### Dependencies

All dependencies should be automatically installed by pip.

###### Installation
- numpy
- scipy
- biopython
- argparse
- pandas
- h5py
- matplotlib (for future use, yet only pylab is used)

### Usage

#### Examples

#### DNA-Example

###### Simulate human DNA reads
```
cd examples/DNA
./Simulator.sh
```
