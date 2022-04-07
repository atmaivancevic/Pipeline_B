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

#### 0a) Append species names to the genome name.
As the names given to genome assemblies are not usually informative, you will want to append species names to the genome names. 

Run [0a_rename_genome.sh](0a_rename_genome.sh). 
Example usage:
```bash
GENOME=<source_genome> SPECIES=<species_name> bash 0a_rename_genome.sh
```

#### 0b) Make each genome a BLAST database and create indexes.
Run [0b_make_database_and_index.sh](0b_make_database_and_index.sh). 
Example usage: 
```bash
bash 0b_make_database_and_index.sh
```

### 1) BLAST TE of interest against all available genomes.

#### 1a) Use TBLASTN with protein sequence queries.
This will identify similar TEs in distantly related species. Output will be nucleotide sequences.
Run [1a_tblastn_and_extract.sbatch](1a_tblastn_and_extract.sbatch).
Example usage:
```bash
DIR=test_genome DATABASE=YarrowiaLipolytica_ASM252v1.fa QUERY=L1_ORFp.fasta RESULTSDIR=results sbatch 1a_tblastn_and_extract.sbatch
```

#### 1b) (Optional) Use BLASTN or LASTZ with nucleotide sequence queries.

Run [1b_lastz_and_extract.sbatch](1b_lastz_and_extract.sbatch).
Example usage:
```bash
GENOMEDIR=test_genome GENOME=YarrowiaLipolytica_ASM252v1.fa QUERYDIR=test_query QUERY=L1_nucl_seqs.fasta RESULTSDIR=results sbatch 1b_lastz_and_extract.sbatch
```

#### 1c) For each genome, combine all identified nucleotide sequences from the previous steps.

Run [1c_combine_hits.sbatch](1c_combine_hits.sbatch).
Example usage:
```bash
SPECIES=YarrowiaLipolytica ELEMENT=L1 LASTZFILE=YarrowiaLipolytica_ASM252v1.fa_L1_nucl_seqs.fasta_lastz.bed TBLASTNFILE=YarrowiaLipolytica_ASM252v1.fa_L1_ORFp.fasta_merged.bed GENOME=YarrowiaLipolytica_ASM252v1.fa RESULTSDIR=results sbatch 1c_combine_hits.sbatch
```

#### 1d) Add header annotations to indicate the genome that each sequence was derived from.

Run [1d_append_name_to_headers.sh](1d_append_name_to_headers.sh).
Example usage:
```bash
SPECIES=YarrowiaLipolytica ELEMENT=L1 RESULTSDIR=results bash 1d_append_name_to_headers.sh
```

Repeat screening in an iterative process (e.g. BLAST-ing the new, larger, query dataset against each genome and then combining the output) until no new hits are found. 

### 2) Perform a reciprocal best hit check.

#### 2a) Use CENSOR to compare hits against known repeat databases (e.g. RepBase or Dfam).
Run [2a_censor_sequences.sbatch](2a_censor_sequences.sbatch).
Example usage:
```bash
INDIR=results FILE=YarrowiaLipolytica_L1_combined.fasta OUTDIR=results/censored sbatch 2a_censor_sequences.sbatch
```

#### 2b) Confirm and extract hits that match the correct TE family.
Run [2b_check_censor_output.sbatch](2b_check_censor_output.sbatch).
Example usage:
```bash
SPECIES=Yarrowia.lipolytica FILE=Yarrowia.lipolytica_L1_combined.fasta GENOME=test_genome/YarrowiaLipolytica_ASM252v1.fa ELEMENT=L1 QUERY=test_query/known_L1_elements_from_repbase.txt CENSORDIR=results/censored sbatch 2b_check_censor_output.sbatch
```

### 3) Cluster all sequences obtained from the iterative alignment screening.

Prior to this step, you will need to combine hits from all genomes into one file. Make sure that sequence headers indicate the species that each TE sequence was derived from.

#### 3a)  All-against-all clustering of nucleotide sequences using VSEARCH.
You can use full-length nucleotide sequences, or nucleotide sequences of the open reading frames only. This clustering step is important as it will reveal likely HTT events which are manifested as clusters of highly similar elements that include elements from multiple species. We have found it best to use sequence divergence cut offs that cluster most closely related sequences (e.g. <20% divergent).

Run [3a_vsearch_cluster_for_nucleotide_seqs.sbatch](3a_vsearch_cluster_for_nucleotide_seqs.sbatch), changing the clustering identity threshold (ID) as required. 
Example usage:
```bash
INDIR=results/allSpeciesCombined FILE=allSpecies_L1.fasta ID=80 PREFIX=c sbatch 3a_vsearch_cluster_for_nucleotide_seqs.sbatch
```

#### 3b) All-against-all clustering of amino acid sequences using USEARCH.
Run [blah](blah).
Example usage:
```bash
bah
```


