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

