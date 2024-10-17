#!/bin/bash
#SBATCH --job-name=example.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=nat64@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#- SET ENVIRONMENT -#
shopt –s expand_aliases
module load trimmomatic

#- SETUP VARIABLES -#
adapters=/home/nat64/Trimmomatic_Monday_Practice/adapters/TruSeq3-PE.fa

input_R1=/home/nat64/Trimmomatic_Monday_Practice/input_files/EC-12_R1_001.fastq.gz
input_R2=/home/nat64/Trimmomatic_Monday_Practice/input_files/EC-12_R2_001.fastq.gz

output_R1_PE=output_R1_pairedend.fq.gz
output_R1_SE=output_R1_singleend.fq.gz

output_R2_PE=output_R2_pairedend.fq.gz
output_R2_PE=output_R2_singleend.fq.gz

#- RUN SLURM -#
trimmomatic PE \
$input_R1 \
$input_R2 \
$output_R1_PE $output_R1_SE \
$output_R2_PE $output_R2_SE \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:10 \
TRAILING:10 \
SLIDINGWINDOW:4:15 \
MINLEN:36

#- UNLOAD MODULE -#
module unload trimmomatic
