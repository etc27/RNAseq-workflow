# Workflow for analyzing differential gene expression from RNA-sequencing data
Emma Tung Corcoran (10/07/2020)

## Introduction
This document covers my basic workflow for processing paired-end RNA sequencing samples and analyzing differential gene expression. I used the [Ruddle HPC cluster at the Yale Center for Research Computing](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/ruddle/) for my HPC environment.

## Setup

### Setting up Miniconda
I followed the directions on the [YCRC Conda Documentation](https://docs.ycrc.yale.edu/clusters-at-yale/guides/conda/) to set up Miniconda as a module on my HPC account. First, I generated a conda environment (name=env_name) containing all of the packages I need for RNA-seq that are not included with the default conda installation.

```
module load miniconda
conda create -n env_name fastqc trim-galore star subread multiqc samtools
```

Then, I am able to load the conda environment containing all of the required packages with the following code.

```
module load miniconda
conda activate env_name
```
### Setting up the Folder Structure
In order to organize all of the files generated from processing the RNA-seq raw data, I utilized the following folder structure adapted from [this RNA-seq workflow](https://github.com/twbattaglia/RNAseq-workflow).

```
── RNAseq_data/
  │   └── annotation/               <- Genome annotation file (.GTF/.GFF)
  │  
  │   └── genome/                   <- Host genome file (.FASTA)
  │  
  │   └── input/                    <- Location of input RNAseq data
  │  
  │   └── output/                   <- Data generated during processing steps
  │       ├── 1_initial_qc/         <- Main alignment files for each sample
  │       ├── 2_trimmed_output/     <-  Log from running STAR alignment step
  │       ├── 3_aligned_sequences/  <- Main alignment files for each sample
  │           ├── aligned_bam/      <-  Alignment files generated from STAR (.BAM)
  │           ├── aligned_logs/     <- Log from running STAR alignment step
  │       ├── 4_final_counts/       <- Summarized gene counts across all samples
  │       ├── 5_multiQC/            <- Overall report of logs for each step
  │  
  │   └── star_index/               <-  Folder to store the indexed genome files from STAR 
  ```
  
