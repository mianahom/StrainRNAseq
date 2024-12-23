#!/bin/bash
#SBATCH --job-name=htseq_count
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=11G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[0-15]%16

hostname
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5
module load parallel/20180122

INDIR=../../results/03_align/alignments
OUTDIR=../../results/05_counts/counts
mkdir -p $OUTDIR

# accession list

SAMPLES=../01_QC/samples.txt
NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)

SAMPLE=$(sed -n ${NUM}p $SAMPLES)

# gtf formatted annotation file
GTF=../../genome/Mus_musculus.GRCm39.113.gtf

# run htseq-count on each sample in parallel
echo $SAMPLE
