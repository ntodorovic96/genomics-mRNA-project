#!/bin/bash
#SBATCH --job-name=example.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=nat64@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#- SET ENVIRONMENT -#
module load bowtie2/2.5.3

#- SETUP VARIABLES -#
ref_genome_index_prefix=/home/nat64/RNAseq_project/ref_genome/bt2-index/ref_gen_index
forward_PE_reads=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_1.trPE.fq.gz
reverse_PE_reads=/home/nat64/RNAseq_project/trimmming/output_files/WTC2_2.trPE.fq.gz

#- RUN SLURM -#
bowtie2 -x $ref_genome_index_prefix \
-1 $forward_PE_reads \
-2 $reverse_PE_reads \
-S WTC2.sam

#- UNLOAD MODULE -#
module unload bowtie2
