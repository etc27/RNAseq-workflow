# Workflow for analyzing differential gene expression from RNA-sequencing data
Emma Tung Corcoran (10/07/2020)

## Introduction
This document covers my basic workflow for processing paired-end RNA sequencing samples and analyzing differential gene expression. I used the [Ruddle HPC cluster at the Yale Center for Research Computing](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/ruddle/) for my HPC environment. Sections of this workflow were adapted from [this RNA-seq tutorial](https://github.com/twbattaglia/RNAseq-workflow).

## 1. Setup

### Set up Miniconda environment
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
### Set up the folder structure
In order to organize all of the files generated from processing the RNA-seq raw data, I utilized the following folder structure adapted from [this RNA-seq workflow](https://github.com/twbattaglia/RNAseq-workflow).

```
── RNAseq_data/
  │   └── annotation/               <- Genome annotation file (.GTF/.GFF)
  │  
  │   └── genome/                   <- Reference genome file (.FASTA)
  │  
  │   └── input/                    <- Location of input RNA-seq data
  │  
  │   └── results/                  <- Data generated during processing steps
  │       ├── 1_initial_qc/         <- Quality check of input files
  │       ├── 2_trimmed_output/     <- Trimmed read files and quality check of trimmed reads
  │       ├── 3_aligned_sequences/  <- Main alignment files for each sample
  │           ├── aligned_bam/      <- Alignment files generated from STAR (.BAM)
  │           ├── aligned_logs/     <- Log from running STAR alignment step
  │       ├── 4_final_counts/       <- Summarized gene counts across all samples
  │       ├── 5_multiQC/            <- Overall report of logs for each step
  │  
  │   └── star_index/               <- Folder to store the indexed genome files from STAR 
  ```
  
### Download the reference genome and annotation
I downloaded the *Arabidopsis thaliana* reference genome (Araport 11) from the [JGI Genome Porta](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Athaliana) to the `genome/` folder. The genome assembly was called `Athaliana_447_TAIR10.fa.gz`
I downloaded the *Arabidopsis thaliana* annotation (Araport 11) from [TAIR](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release) to the `annotation/` folder. The annotation was called `Araport11_GFF3_genes_transposons.201606.gtf`

### Download raw sequencing data
In order to access the sequencing data from the RNA-seq experiments, I followed the directions on the [Ruddle documentation](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/ruddle/#access-sequencing-data). Briefly, the Yale Center for Genome Analysis sent me a url that looks like this: 
`http://fcb.ycga.yale.edu:3010/randomstring/sample_dir_001` and I used the ycgaFastq tool to make a soft link to the data with the following command.
```
/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq  fcb.ycga.yale.edu:3010/randomstring/sample_dir_001
```
My raw sequencing data contains paired-end Illumina sequencing reads that were rRNA-depleted prior to sequencing. Note that the fastq files corresponding to pair1 and pair2 for one sample are labeled `sample_R1_001.fastq.gz` and `sample_R2_001.fastq.gz`

## 2. Analyze sequence quality with FastQC

### Description
[FastQC: A quality control tool for high throughput sequence data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis."

### Command
Note: it is not necessary to unzip the fastq file for any of the following processes.
```
# Run FastQC
# -o: output directory
# --noextract: do not uncompress the output file after creating it

fastqc -o results/1_initial_qc/ --noextract input/sample.fastq.gz
```

### Output
```
── results/1_initial_qc/
    └──  sample_fastqc.html   <- HTML file of FastQC quality analysis figures
    └──  sample_fastqc.zip    <- FastQC report data
```

## 3. Perform quality control with Trim Galore

### Description
[Trim Galore: A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
Trim Galore performs adapter trimming (the default is the first 13 bp of Illumina standard adapters ('AGATCGGAAGAGC'), but it is able to autodetect the adapter sequence). It also removes sequences that become too short during the trimming process. With the `--paired` option, Trim Galore removes both reads in a pair if at least one of the two sequences becomes shorter than the threshold. Additionally, it can run FastQC on the output files to assess quality once trimming has been completed. I kept the default options for quality and length.

### Command
```
# Run Trim Galore! (input is both paired end sequencing files for a sample)
#--paired: remove both reads in a pair if at least one of the two sequences becomes shorter than the threshold
#--fastqc: run FastQC in the default mode on the FastQ file once trimming is complete
#--output_dir: output directory

trim_galore --paired --fastqc --output_dir results/2_trimmed_output/ input/sample_R1_001.fastq.gz input/sample_R2_001.fastq.gz
```

### Output
```
── results/2_trimmed_output/
     └──  sample_R1_001_val_1.fq.gz                   <- Compressed trimmed sequencing file (for read1)
     └──  sample_R1_001_val_1_fastqc.html                    <- HTML file of FastQC quality analysis figures (for read1)
     └──  sample_R1_001_val_1_fastqc.zip                     <- FastQC report data (for read1)
     └──  sample_R1_001.fastq.gz_trimming_report.txt  <- Cutadapt trimming report (for read1)
     └──  sample_R2_001_val_2.fq.gz                   <- Compressed trimmed sequencing file (for read2)
     └──  sample_R2_001_val_2_fastqc.html                    <- HTML file of FastQC quality analysis figures (for read2)
     └──  sample_R2_001_val_2_fastqc.zip                     <- FastQC report data (for read2)
     └──  sample_R2_001.fastq.gz_trimming_report.txt  <- Cutadapt trimming report (for read2)
```

## 4. Align reads to genome with STAR

### Description
[STAR: ultrafast universal RNA-seq aligner](https://pubmed.ncbi.nlm.nih.gov/23104886/)
STAR (Spliced Transcripts Alignment to a Reference) is an ultrafast alignment software for aligned RNA-seq data to genomes.

### Generating index
Before running STAR, it is necessary to generate an index of your reference genome. I stored this index in the `star_index` folder. It is only necessary to run this step once.
```
# Generate genome index for STAR
#--runMode: directs STAR to run genome indices generation job
#--genomeDir: specifies path to the directory where the genome indices are stored
#--genomeFastaFiles: specifies one or more FASTA files with the genome reference sequences
#--sjdbGTFfile: specifies the path to the file with annotated transcripts in the standard GTF format
#--genomeSAindexNbases: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches
#--runThreadN: number of threads

STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome/* --sjdbGTFfile annotation/* --genomeSAindexNbases 12 --runThreadN 4
```

### Command
```
# Run STAR
#--genomeDir: specifies path to the directory where the genome indices are stored
#--readFilesCommand: indicates to uncompress .gz files
#--readFilesIn: path of the files with the sequences to be mapped (for paired-end reads, read1 and read2 files have to be supplied)
#--outFilterMismatchNmax: maximum number of mismatches per pair
#runThreadN: number of threads
#--outSAMtype: output format (BAM SortedByCoordinate = output sorted by coordinate)
#--quantMode: GeneCounts specifies for STAR to count the number of reads per gene while mapping
#--outFileNamePrefix: output files name prefix

STAR --genomeDir star_index --readFilesCommand zcat \
--readFilesIn results/2_trimmed_output/sample_R1_001_val_1.fq.gz results/2_trimmed_output/sample_val_2.fq.gz \
--outFilterMismatchNmax 2 --runThreadN 4 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
--outFileNamePrefix results/4_aligned_sequences/sample
```

### Output
```
── results/3_aligned_sequences/
    └── aligned_bam/sampleAligned.sortedByCoord.out.bam   <- Sorted BAM alignment
    └── aligned_logs/sampleLog.final.out                  <- Summary mapping statistics after mapping job is complete
    └── aligned_logs/sampleLog.out                        <- Main log file with a lot of detailed information about the run
    └── aligned_logs/sampleLog.progress.out               <- Reports job progress statistics taken in 1-minute intervals
    └── aligned_logs/sampleReadsPerGene.out.tab           <- Output containing counts per gene in tab-delimited format
    └── aligned_logs/sampleSJ.out.tab                     <- Contains high confidence collapsed splice junctions in tab-delimited format
```

### Move output files to correct subdirectories
```
# Move the BAM file into the correct folder
mv -v results/3_aligned_sequences/sampleAligned.sortedByCoord.out.bam results/3_aligned_sequences/aligned_bam/

# Move the logs into the correct folder
mv -v results/3_aligned_sequences/*.out results/3_aligned_sequences/aligned_logs/
mv -v results/3_aligned_sequences/*.out.tab results/3_aligned_sequences/aligned_logs/
```

## 5. Summarize gene counts with featureCounts

### Description
[featureCounts: an efficient general purpose program for assigning sequence reads to genomic features](https://pubmed.ncbi.nlm.nih.gov/24227677/)
featureCounts is a read summarization program that can count reads from RNA sequencing experiments (it is also capable of counting reads from DNA sequencing experiments).

### Command
```
# Change directory into the aligned .BAM folder
cd results/3_aligned_sequences/aligned_bam

# Store list of files as a variable
dirlist=$(ls -t ./*.bam | tr '\n' ' ')
echo $dirlist

# Run featureCounts on all of the samples
#-a: path to annotation
#-o: path to output results
#-g: attribute type (i.e. gene_id or gene_name)
#-T: number of threads
featureCounts -a ../../annotation/* -o ../../results/4_final_counts/final_counts.txt -g 'gene_name' -T 4 $dirlist
```

### Output
```
── results/4_final_counts/
    └── final_counts.txt                <- Final gene counts across all samples
    └── final_counts.txt.summary        <- Summary of gene counts
```

## 6. Generate analysis report with MultiQC

### Description
[MultiQC: summarize analysis results for multiple tools and samples in a single report](https://pubmed.ncbi.nlm.nih.gov/27312411/)
MultiQC is a tool to create a single report visualizing output from multiple tools across many samples, including the output of FastQC, Trim_Galore, STAR, and featureCounts.

### Command
```
# Run multiQC
#--outdir: output directory
multiqc results --outdir results/5_multiQC
```

### Output
```
── results/5_multiQC/
    └── multiqc_report.html     <- Figures representing the logs from each step
    └── multiqc_data/           <- Folder of data that multiqc found from various log files
```
