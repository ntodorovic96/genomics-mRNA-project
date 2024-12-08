#!/bin/bash
#SBATCH --job-name=example.SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=nat64@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=24:00:00
#SBATCH --mem=4G

#- SET ENVIRONMENT -#
source activate htseq

#- SETUP VARIABLES -#
htseq=/home/nat64/.conda/envs/htseq/bin/htseq-count
bam_file=/home/nat64/RNAseq_project/seq_alignment/WTC2.srt.bam
gtf_file=/home/nat64/RNAseq_project/seq_alignment/GCF_000182965.3_ASM18296v3_genomic.gtf
output_file=/home/nat64/RNAseq_project/htseq/WTC2_htseq-count.txt

#- RUN SLURM -#
${htseq} -f bam \
-r pos \
-s no \
-t exon \
-i gene_id \
$bam_file \
$gtf_file \
> $output_file
