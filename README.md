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
- Some level of familiarity with computers and queuing systems (e.g. SLURM)

#### Test dataset
A test genome (fungus *Yarrowia lipolytica*) has been placed in [test_genome](test_genome), along with a set of L1 repeats as a [test_query](test_query). We recommend trying out the workflow below on this genome first. 

## Usage

### 1) Extraction of L1 and BovB retrotransposons from genome data

#### 1a) Download or acquire genomes of interest
755 genomes were downloaded from public databases (UCSC and NCBI); four more were acquired from collaborators. All genomes were downloaded using ```wget```, as recommended by NCBI (https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/). See Supplementary Table 1 for the source and assembly version of each genome used.

#### 1b) Append species name to each genome fasta file
This will save you a lot of headaches later on. Most genomes are given useless names like GCF_000002525.2_ASM252v1_genomic.fna.
Use [renameGenome.q](scripts/renameGenome.q). 


