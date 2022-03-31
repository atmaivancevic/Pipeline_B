#!/bin/bash

# Index all genomes and make BLAST database
# bash 0_make_database_and_index.sh

# Load modules if you are working on a SLURM enivronment, e.g. 
module load BLAST+
module load SAMtools

# Make each genome a BLAST database

for i in genomes/*.fa; do makeblastdb -in $i -parse_seqids -dbtype nucl; done

# Create a fai index for each genome

for i in genomes/*.fa; do samtools faidx $i; done
