#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################
# this time we'll use qualimap

# load software--------------------------------------------------------------------------
module load qualimap/2.3
module load parallel/20180122

# input, output directories--------------------------------------------------------------

INDIR=../../results/03_align/alignments
OUTDIR=../../results/04_alignQC/qualimap_reports
mkdir -p $OUTDIR

# gtf annotation is required here
GTF=../../Mus_musculus.GRCm39.113.gtf

# accession list
SAMPLES=../01_QC/samples.txt

# run qualimap in parallel
cat $SAMPLES | \
parallel -j 5 \
    qualimap \
        rnaseq \
        -bam $INDIR/{}.bam \
        -gtf $GTF \
        -outdir $OUTDIR/{} \
        --java-mem-size=10G  
