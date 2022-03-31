#!/bin/bash

# Index all genomes and make BLAST database
# bash 0_make_database_and_index.sh

# Make each genome a BLAST database
module load BLAST+/2.2.31-foss-2015b-Python-2.7.11
for i in genomes/*.fna; do makeblastdb -in $i -parse_seqids -dbtype nucl; done

# Create a fai index for each genome
module load SAMtools/1.3.1-foss-2016b
for i in genomes/*.fna; do samtools faidx $i; done
