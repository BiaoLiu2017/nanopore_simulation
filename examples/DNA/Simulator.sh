#!/bin/bash

# if [ ! -f ../Homo_sapiens.GRCh38.dna.primary_assembly.fa ]; then
#     wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#     gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ../Homo_sapiens.GRCh38.dna.primary_assembly.fa
# fi
# if [ -d Run-Output ]; then
#     rm -r Run-Output
# fi
# conda activate Simulator
python /home/liubiao/Chiron/Simulator/nanopore_simulation/bin/simulatION.py \
simulate \
-r /home/liubiao/Chiron/Simulator/nanopore_simulation/examples/DNA/test.fasta \
-c /home/liubiao/Chiron/Simulator/nanopore_simulation/examples/DNA/config.txt \
-m /home/liubiao/Chiron/Simulator/nanopore_simulation/examples/DNA/BGI_4mer.txt \
-n 100 \
-o /home/liubiao/Chiron/Simulator/nanopore_simulation/Run-Output \
--error_rate 0.0 \
--reverse False \
--signal_repeat True \
--uniform_noise True \
--event_repeat True \
--low_pass_filter True \
--gaussian_noise True \
--correction True
               
               
