#!/bin/bash
#SBATCH --job-name=example.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=nat64@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#- SET ENVIRONMENT -#
shopt -s expand_aliases
module load trimmomatic

#- SETUP VARIABLES -#
adapters=/home/nat64/RNAseq_project/input_files/adapters/TruSeq3-PE.fa

input_R1=/home/nat64/RNAseq_project/input_files/PE_reads/WTC2_1.fq.gz
input_R2=/home/nat64/RNAseq_project/input_files/PE_reads/WTC2_2.fq.gz

output_R1_PE=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_1.trPE.fq.gz
output_R1_SE=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_1.trSE.fq.gz

output_R2_PE=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_2.trPE.fq.gz
output_R2_SE=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_2.trSE.fq.gz

#- RUN SLURM -#
trimmomatic PE \
$input_R1 \
$input_R2 \
$output_R1_PE $output_R1_SE \
$output_R2_PE $output_R2_SE \
ILLUMINACLIP:$adapters:2:30:10 \
HEADCROP:15 \
TRAILING:20 \
SLIDINGWINDOW:4:15 \
MINLEN:75

#- UNLOAD MODULE -#
module unload trimmomatic
