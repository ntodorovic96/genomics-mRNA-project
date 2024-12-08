#!/bin/bash
#SBATCH --job-name=example.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=nat64@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#- SET ENVIRONMENT -#
shopt -s expand_aliases
module load fastqc

#- SETUP VARIABLES -#
output_folder=/home/nat64/RNAseq_project/trimmming/fastqc_ouputs
fastq_file_1=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_1.trPE.fq.gz
fastq_file_2=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_2.trPE.fq.gz

#- RUN SLURM -#
fastqc -o $output_folder $fastq_file_1 $fastq_file_2

#- UNLOAD MODULE -#
module unload fastqc
