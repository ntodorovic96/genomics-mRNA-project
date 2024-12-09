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

To resolve these issues, Trimmimotatic was used to remove any adapter content, remove the first 15 bases, eliminate trailing bases with a low quality score, remove potential reads with low average quality socres and remove any reads that were less than 75 basepairs, in that order (main/scripts/trimming.SBATCH). The input fastq files used were WTC2_1.fq.gz and WTC2_2.fq.gz, and the files generated were WTC2_1.trPE.fq.gz and WTC2_2.trPE.fq.gz. These trimmed files were then ran through FastQC to analyze the results of the trimming, generating two html files WTC2_1.trPE_fastqc.html and WTC2_2.trPE_fastqc.html (main/scripts/fastqc2.SBATCH). Notably, the per-base sequence quality issues were resolved due to the heacrop, and adapter content was completely removed. Even though FastQC raised concerns about sequence duplication levels in both the precleaned and cleaned data, we didn't remove duplicate reads because duplicates indicate potentially useful information for expression levels in RNAseq. At the end of the cleaning, were 19,434,301 reads per library and 95% of the reads were retained (main/spreadsheets/trimming_data.csv).

# Reference Genome Mapping
To map the RNAseq reads, we utilized a _C. Albicans_ reference genome downloaded from NCBI (GCF_000182965.3). The fasta file used from this genome assembly was GCF_000182965.3_ASM18296v3_genomic.fna and the annotation file used was GCF_000182965.3_ASM18296v3_genomic.gtf.

Using Bowtie-2, the fasta reference genome file was indexed to generate an indexed map (main/additional-commands/bt2-build). Once the indexed map was generated, Bowtie-2 was used to map the trimmed reads against the reference genome (main/scripts/bt2.SBATCH). The input files into the script were the trimmed paired-end reads WTC2_1.trPE.fq.gz and WTC2_2.trPE.fq.gz, and the indexed reference genome ref_gen_index. The output of this alignment was the sequence alignment map (SAM) file WTC2.sam. Out of the 19,434,301 reads, all of them were detected and 98.07% were aligned (main/spreadsheets/bowtie_alignments.csv). Additionally, 89.12% of the reads were aligned exactly 1 time, while 5.91% were aligned more than once. 

Finally, the SAM file was converted into a binary alignment map (BAM) file WTC2.bam (main/additional-commands/bam_conversion) and sorted as well as indexed (main/additional-commands/bam-srt). 

# Read Counts (HTSeq)
After mapping the reads to a BAM file with a reference genome, the python program HTSeq was used to count the reads. To begin, a Conda virtual environment (VE) was set up with all of the dependencies and necessary packages to run HTSeq (main/additional-commands/create-conda). Subsequently, HTSeq was run in the conda VE, using the WTC2.srt.bam sorted BAM file, and GCF_000182965.3_ASM18296v3_genomic.gtf GTF file as inputs (main/scripts/htseq.SBATCH). The generated output file was the text file WTC2_htseq-count.txt.

# Differential Expression Analysis
All of the 6 HTseq read count files from each biological replicate in thiamine-present and thiamine absent conditions (WTA2_htseqCount, WTB2_htseqCount, WTC2_htseqCount, WTA1_htseqCount, WTB1_htseqCount, WTC1_htseqCount) were then inputted into a DESeq R script (main/scripts/DESeq_script.R).

DESeq made a dataset from all the HTseq count data and filtered for the relevant genes with more than 10 reads across samples (main/spreadsheets/calb_TH-_v_TH+_allgenes.csv). The reference condition to base the fold change in the differential expression analysis was set to the thiamine-present condition, and each experimental group and replicate was separated based on two principle components into a principle component analysis (PCA) plot. Finally, the data was filtered down further to a set of differentially expressed genes (DEGs) with at least 2-fold increase in expression and an adjusted p-value of less than 0.05 (main/spreadsheets/signif_TH-vTH+_wnames.csv). This subset of DEGs were graphed on a volcano plot.

# Results of DEG Analysis
The pre-cleaning for genes that had at least ten reads across samples resulted in the identification of 6,072 potential gene candidates.
![TH-vTH+_pcaplot](https://github.com/user-attachments/assets/a82174d6-d051-4e7a-8a67-f4b7109a18d5)
Using these candidate genes and their expression levels, each biological replicate and treatment group was separated based on two principle components. Principle Component 1 (PC1) accounted for 88% of the observed variance, while Principle Component 2 (PC2) accounted for 9% of the observed variance, leaving only 3% of the variance unaccounted for on the PCA plot. 

The treatment conditions were very well separated based on PC2, with all of the thiamine-absent samples scoring very low at a -10, while all of the thiamine-present samples all scored highly around a +10, indicating some differential expression of genes based on the conditions of thiamine in the environment. On the other hand, PC1 separated all of the biological replicates, with two scoring high in each experimental group and one scoring low. 

<img width="416" alt="Screenshot 2024-12-08 180736" src="https://github.com/user-attachments/assets/ca572482-924b-4525-879b-b48b2d71f768">

Following the filtering down of candidates based on a combination of p-score and fold change criteria, 13 DEGs were identified in thiamine-absent conditions. Interestingly, all 13 of these DEGs were upregulated in thiamine-absent conditions, indicated by the the positive log-fold change values. 

# Physiological Function of DEGs and Gene Ontology Analysis
To identify the physiological role of the identified DEGs, the locus tags of each candidate were pulled from the R output table and used to parse the _C. Albicans_ reference GTF file downloaded earlier (GCF_000182965.3_ASM18296v3_genomic.gtf) for significant gene annotation information, including the db_xref qualifier. These genes were then researched to identify physiological function and involvement in cellular processes on https://www.candidagenome.org/.

<img width="629" alt="Screenshot 2024-12-08 223836" src="https://github.com/user-attachments/assets/0abd7945-37d9-4c60-9574-8f526562e75d">

To better interpret the DEGs in the context of one another and identify the biological processes that these genes are overarchingly involved in, a Gene Ontology (GO) Enrichment Analysis was then conducted. Results were only displayed for a False Discovery Rate (FDR) of P < 0.05 using a Fischer's Exact Test.

<img width="392" alt="Screenshot 2024-12-08 224511" src="https://github.com/user-attachments/assets/15b35422-ad3b-45a6-a291-60bfbe20a341">


