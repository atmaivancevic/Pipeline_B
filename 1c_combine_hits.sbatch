#!/bin/bash

# Example usage:
# SPECIES=human ELEMENT=L1 LASTZFILE=hg38_lastz_output.bed TBLASTNFILE=hg38_tblastn_output.bed GENOME=hg38.fa RESULTSDIR=. sbatch 1c_combine_hits.sbatch

#SBATCH -p batch
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=12:00:00 
#SBATCH --mem=32GB    

# Load necessary modules
module load BEDTools/2.25.0-foss-2015b

# Combine bed files
cat $TBLASTNFILE $LASTZFILE > "$SPECIES"_"$ELEMENT"_together.tmp

# Sort and merge hits
bedtools sort -i "$SPECIES"_"$ELEMENT"_together.tmp > "$SPECIES"_"$ELEMENT"_together.sorted.tmp
bedtools merge -s -i "$SPECIES"_"$ELEMENT"_together.sorted.tmp -c 4 -o collapse > "$SPECIES"_"$ELEMENT"_together.merged.tmp

# Reformat
awk '{print $1 "\t" $2 "\t" $3 "\t" "BovB" "\t" "1" "\t" $4}' "$SPECIES"_"$ELEMENT"_together.merged.tmp > $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed

# Extract nucleotide sequences from bed file
bedtools getfasta -s -fi $GENOME -bed $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed -fo $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.fasta

# Remove tmp files
rm *tmp
