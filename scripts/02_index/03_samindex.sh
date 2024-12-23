#!/bin/bash
#SBATCH --job-name=samindex
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
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
# Download genome and annotation from ENSEMBL
#################################################################

# load software
module load samtools/1.16.1

cd ../../genome
# generate simple samtools fai indexes 
samtools faidx Mus_musculus.GRCm39.dna.toplevel.fa
