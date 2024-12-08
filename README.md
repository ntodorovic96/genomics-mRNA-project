# Candida Albicans RNA-seq project
This project is part of a genomics and bioinformatics class at Georgetown University with the objective of analyzing RNA-seq data from _Candida Albicans_.

# Project Description
The question we are looking to answer is why C. Albicans, a normally commensal yeast, can become pathogenic in certain environmental conditions and causes various human pathologies, such as urinary tract, genital, mucosal and blood infections.

One hypothesis for why this transition may be occurring is that they only become pathogenic in nutrient-deficient environments. To test this hypothesis, we have grown C. Albicans in Thiamine-present and Thiamine-absent media and used RNAseq to analyze their transcriptomic profiles in these conditions. 

# Data Collection and Raw Data
Original FastQ file sequencing data was generated by Dr. Ronda Rolfes and generously allowed us to use it. 
The data collected was paired-end sequencing data, and three biological replicates were made per treatment group. 

From the experimental workflow, 12 raw data files were generated:
WTA1_1.fq.gz, WTA1_2.fq.gz
WTB1_1.fq.gz, WTB1_2.fq.gz
WTC1_1.fq.gz, WTC1_2.fq.gz

WTA2_1.fq.gz, WTA2_2.fq.gz
WTB2_1.fq.gz, WTB2_2.fq.gz
WTC2_1.fq.gz, WTC2_2.fq.gz

The naming conventions are as follows: A/B/C correspond to biological replicates, A1/B1/C1 are from the Thiamine-present group, A2/B2/C2 are from the Thiamine-absent group, and _1,_2 are matched read pairs. 

# Trimming and Data Cleaning
The FastQ data was processed by initially running it through FastQC to guide the subsequent trimming strategy (main/scripts/fastqc.SBATCH). The input fastq files used were WTC2_1.fq.gz and WTC2_2.fq.gz, and the files generated were WTC2_1_fastq.html and WTC2_2_fastq.html. The main areas of concern with regards to sequencing quality were the per-base sequence quality and the per-tile sequence quality. When looking at the per-base sequence quality, the fluctuations in sequence content mainly came from the first 10-15 bases. Additionally, there was a small amount of adapter content, even though FastQC didn't raise any serious concerns. 

To resolve these issues, Trimmimotatic was used to remove any adapter content, remove the first 15 bases, eliminate trailing bases with a low quality score, remove potential reads with low average quality socres and remove any reads that were less than 75 basepairs, in that order (main/scripts/trimming.SBATCH). The input fastq files used were WTC2_1.fq.gz and WTC2_2.fq.gz, and the files generated were WTC2_1.trPE.fq.gz and WTC2_2.trPE.fq.gz. These trimmed files were then ran through FastQC to analyze the results of the trimming, generating two html files WTC2_1.trPE_fastqc.html and WTC2_2.trPE_fastqc.html. Notably, the per-base sequence quality issues were resolved due to the heacrop, and adapter content was completely removed. Even though FastQC raised concerns about sequence duplication levels in both the precleaned and cleaned data, we didn't remove duplicate reads because duplicates indicate potentially useful information for expression levels in RNAseq. At the end of the cleaning, were 19,434,301 reads per library and 95% of the reads were retained (main/spreadsheets/trimming_data.csv).

# Reference Genome Mapping
To map the RNAseq reads, we utilized a _C. Albicans_ reference genome downloaded from NCBI (GCF_000182965.3). The fasta file used from this genome assembly was GCF_000182965.3_ASM18296v3_genomic.fna and the annotation file used was GCF_000182965.3_ASM18296v3_genomic.gtf. Using Bowtie-2, the fasta reference genome file was indexed to generate an indexed map (no script used here, simply bowtie2-build _reference-genome_ _index_). Once the indexed map was generated, Bowtie-2 was used to map the trimmed reads against the reference genome (main/scripts/bt2.SBATCH). The input files into the script were the trimmed paired-end reads WTC2_1.trPE.fq.gz and WTC2_2.trPE.fq.gz, and the indexed reference genome ref_gen_index. The output of this alignment was the sequence alignment map (SAM) file WTC2.sam. 
