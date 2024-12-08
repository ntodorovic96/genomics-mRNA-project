
########################## Installations ###########################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")
BiocManager::install("GOplot")
BiocManager::install("mygene")
install.packages("RColorBrewer")
install.packages("ggolot2")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("viridis")
install.packages("ggthemes")
install.packages("VennDiagram")
BiocManager::install("genefilter")
BiocManager::install("ggrepel")


########################## Load Libraries ###########################
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(GOplot) 
library("apeglm")
library(mygene)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(DESeq2)
library(VennDiagram)
library(genefilter)
library(ggrepel)

########################## Input HTSeq data files ###########################
#Choose directory with htseq-count data
directory<-("/Users/paa9/Desktop/htseq_counts/")
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- list.files(directory)
head(sampleFiles)
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
#sampleNames <- substr(sampleNames, 1, nchar(sampleNames)-4) # this keeps only the treatment plus the replicate
head(sampleNames)
sampleConditions <- substr(sampleFiles, 1, 3)#to get conditions I'm pulling the first letter, which is either N (NTH) or T (TH_)
head(sampleConditions)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions) ## condition is either N for (NTH) or T for TH
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)
View(sampleTable)

########################## Make the DESeq dataset from this HTSeq count data ###########################

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds
View(dds)

########################## Pre-filtering ###########################
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. 
#They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

########################## Set the reference condition ###########################
# Set the condition for what to compare to. 
## (Default is first condition alphabetically)
dds$condition <- relevel(dds$condition, ref = "TH_") # setting reference condition as TH_ which is the THI+ group
head(dds$condition)


########################## Quality control then PCA Visualization ##########################
cds <- estimateSizeFactors(dds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds, blind=TRUE)
theme_set(theme_bw())
#meanSdPlot(assay(vsd))
#plotDispEsts(cds)
nudge <- position_nudge(y = 10)
PCA_data <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE,ntop=500)
View(PCA_data)
percentVar <- round(100 * attr(PCA_data, "percentVar"))
nudge <- position_nudge(y = 4)

pca_plot <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = condition, position_nudge(y=10))) +
  geom_point(size =3, position = position_jitter(w=0.05, h=0.05)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("dark blue", "dark red"),labels=c("THI+","THI-"))+
  labs(color = "Treatment")
# geom_text_repel()
# By default plotPCA() uses the top 500 most variable genes. 
# You can change this by adding the ntop= argument and specifying how many of the genes you want the function to consider.
setwd("/Users/paa9/Desktop/DESeq2_analysis/")
ggsave("TH-vTH+_pcaplot.png",plot=pca_plot,dpi=600,units='in',width=6,height=5)
pca_plot

########################## Run DEG Analysis ###########################
data <- DESeq(dds)

# DESeq2 uses a negative binomial distribution to model RNA-seq counts since they exhibit overdispersion
# (variance > mean). 
# It is common to shrink the log fold change estimate for better visualization and ranking of genes.
#Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(data)
res_LFC <- lfcShrink(data, coef="condition_NTH_vs_TH_", type="apeglm") ## plug in the output from resultsNames(dds) line as the coef
res_LFC
head(res_LFC)
# Order the table by smallest p value
resOrdered <- res_LFC[order(res_LFC$padj),]
summary(resOrdered)
head(resOrdered)
View(resOrdered)

# Write out a table of all genes
setwd("/Users/paa9/Desktop/DESeq2_analysis/")
write.csv(resOrdered, 
          file="./calb_TH-_v_TH+_allgenes.csv", row.names = T)

# Make a volcano plot
## visualizer each results set with a volcano plot
VP <- EnhancedVolcano(res_LFC,
                      lab = rownames(res_LFC),
                      x = 'log2FoldChange',
                      y = 'padj',
                      ylab = "-Log10(p-adjusted)",
                      selectLab = NA,
                      #drawConnectors = TRUE,
                      xlim = c(-2.5, 2.5),
                      ylim = c(0,25),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 2.0,
                      labSize = 5.0)
VP

ggsave("volcano_TH-vTH+.png",plot=VP,dpi=600,units='in',width=8,height=6)

sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)
View(sig_res)
length(unique(sig_res$gene)) #13 


# adding gene names to this deg
genes <- queryMany(sig_res$gene, fields=c("name"), species="237561")
genes
colnames(genes)<-c("gene","id","score","gene name")
genes <- as.data.frame(subset(genes, select = c("gene", "gene name")))
View(genes)
length(unique(genes$gene))
length(genes$gene)
genes <- (genes[!duplicated(genes), ]) ##getting rid of duplicate gene name rows for Trnad-guc
sig_res.gene <- as.data.frame(merge(genes,sig_res,by="gene"))
View(sig_res.gene) 
length(unique(sig_res.gene$gene))
length(sig_res.gene$gene)

# Write out a table of these significant differentially expressed genes
setwd("/Users/paa9/Desktop/DESeq2_analysis/")
write.csv(sig_res.gene, 
          file="./signif_TH-vTH+.csv", row.names = F)

