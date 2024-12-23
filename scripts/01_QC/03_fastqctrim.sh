#!/bin/bash
#SBATCH --job-name=fastqc_trim
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastqc/0.12.1
module load parallel/20180122
module load MultiQC/1.15

INDIR=../../results/01_QC/trim_reads
REPORTDIR=../../results/01_QC/trimmed/fastqc
mkdir -p $REPORTDIR
MULTIQC=../../results/01_QC/trimmed/multiqc
mkdir -p $MULTIQC


# run fastqc in parallel, 6 samples at a time
ls $INDIR/*_trim_{1..2}.fastq.gz | parallel -j 6 \
    fastqc --outdir $REPORTDIR {}



# run on fastqc output
multiqc -f -o $MULTIQC $REPORTDIR
