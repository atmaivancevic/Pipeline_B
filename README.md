# Pipeline B

This workflow is used for detecting potential horizontal transfer events, starting from a set of curated repeat consensus sequences from available sources (e.g. RepBase or Dfam) for the TE of interest. 

It is a simplified version of the code used to infer horizontal transfer events involving L1 and BovB retrotransposons in eukaryotes, as described here: https://github.com/AdelaideBioinfo/horizontalTransfer

#### Requirements:
- BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- LASTZ (https://www.bx.psu.edu/~rsharris/lastz/)
- SAMtools (http://www.htslib.org/)
- BEDtools (http://bedtools.readthedocs.io/en/latest/)
- CENSOR, which requires wu-blast and bioperl (https://girinst.org/downloads/software/censor/)
- VSEARCH (https://github.com/torognes/vsearch)
- USEARCH (https://www.drive5.com/usearch/)

#### Prerequisites
- Some level of familiarity with computers and queuing systems (e.g. SLURM).

#### Test dataset
A test genome (fungus *Yarrowia lipolytica*) has been placed in [test_genome](test_genome), along with a set of L1 repeats as a [test_query](test_query). We recommend trying out the workflow using these files first.

## Usage

### 0) Download and prepare the genomes you want to screen.

#### 0a) Append species names to the genome name
As the names given to genome assemblies are not usually informative, you will want to append species names to the genome names. 

Run [0a_rename_genome.sh](0a_rename_genome.sh). 
Example usage:
```bash
GENOME=<source_genome> SPECIES=<species_name> bash 0a_rename_genome.sh
```

#### 0b) Make each genome a BLAST database and create indexes
Run [0b_make_database_and_index.sh](0b_make_database_and_index.sh). 
Example usage: 
```bash
bash 0b_make_database_and_index.sh
```

### 1) BLAST TE of interest against all available genomes.

#### 1a) Use TBLASTN with protein sequence queries
This will identify similar TEs in distantly related species. Output will be nucleotide sequences.
Run [1a_tblastn_and_extract.sbatch](1a_tblastn_and_extract.sbatch).
Example usage:
```bash
DIR=test_genome DATABASE=YarrowiaLipolytica_ASM252v1.fa QUERY=L1_ORFp.fasta RESULTSDIR=results sbatch 1a_tblastn_and_extract.sbatch
```

#### 1b) (Optional) Use BLASTN or LASTZ with nucleotide sequence queries

Run [1b_lastz_and_extract.sbatch](1b_lastz_and_extract.sbatch).
Example usage:
```bash
GENOMEDIR=test_genome GENOME=YarrowiaLipolytica_ASM252v1.fa QUERYDIR=test_query QUERY=L1_nucl_seqs.fasta RESULTSDIR=results sbatch 1b_lastz_and_extract.sbatch
```

#### 1c) For each genome, combine all identified nucleotide sequences from the previous steps

Run [1c_combine_hits.sbatch](1c_combine_hits.sbatch).
Example usage:
```bash
SPECIES=YarrowiaLipolytica ELEMENT=L1 LASTZFILE=YarrowiaLipolytica_ASM252v1.fa_L1_nucl_seqs.fasta_lastz.bed TBLASTNFILE=YarrowiaLipolytica_ASM252v1.fa_L1_ORFp.fasta_merged.bed GENOME=YarrowiaLipolytica_ASM252v1.fa RESULTSDIR=results sbatch 1c_combine_hits.sbatch
```

#### 1d) Add header annotations to indicate the genome that each sequence was derived from

Run [1d_append_name_to_headers.sh](1d_append_name_to_headers.sh).
Example usage:
```bash
SPECIES=YarrowiaLipolytica ELEMENT=L1 RESULTSDIR=results bash 1d_append_name_to_headers.sh
```



