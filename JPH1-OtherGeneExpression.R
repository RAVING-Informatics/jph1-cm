#Generate a plot showing epxression of RYR1, CACNA1S, JPH2, JPH3, MTM1, DNM2, BIN1, STAC3, ORAI1, STIM1, CAV3 and TTN 
#in patient F3-II:1 compared other NMD patients and controls.

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofWesternAustralia/PERKINS/RNAseq/muscle_counts/")

#load packages
library(DESeq2)
library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(OUTRIDER)

## using normalised count data from OUTRIDER data set
ods <- readRDS(file = "/Users/00104561/Library/CloudStorage/OneDrive-TheUniversityofWesternAustralia/Research/Projects/RNAseq/DROP/Muscle_May2023/processed_results/aberrant_expression/v38/outrider/MUSCLE/ods.Rds")
counts_norm<-as.data.frame(counts(ods, normalized = TRUE))
counts_norm$geneID<-row.names(counts_norm)

#lookup the correct Ensembl version using IDs
rownames(counts_norm)[grep("ENSG00000154118", rownames(counts_norm))]

#create a gene list for RYR1, CACNA1S, JPH2, JPH3, MTM1, DNM2, BIN1, STAC3, ORAI1, STIM1, CAV3 and TTN
genes<-c("ENSG00000104369.5","ENSG00000196218.14", "ENSG00000149596.7", "ENSG00000171100.16", "ENSG00000079805.19", "ENSG00000136717.15", "ENSG00000185482.8", "ENSG00000276045.3", "ENSG00000167323.12", "ENSG00000182533.7", "ENSG00000155657.29", "ENSG00000081248.12")
genes

#reorder counts for plotting 
reord_counts<- counts_norm[genes,] %>%
  gather("sample_ID", "counts", 0:129)

##Get a list of hgnc_symbols based on the ensembl gene ids in count data
#use biomart to retrieve a list of gene names using the ensembl gene ID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- reord_counts$geneID
genes_sub <- sub("\\.\\d+$", "", genes) #remove the version from the gene_ID
reord_counts$genes_sub<-genes_sub #add ensembl id without version number
G_list <- getBM(filters="ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=genes_sub,mart= mart)

#merge the hgnc gene names with the original results table
counts_hgnc <- unique(merge(reord_counts,G_list,by.x="genes_sub",by.y="ensembl_gene_id", all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE))
head(counts_hgnc) #check

##plots for each gene of interest, highlighting the patient (D20-0017)
highlight <- counts_hgnc %>%
  filter(sample_ID=="D20-0017")
other <- counts_hgnc %>%
  filter(sample_ID!="D20-0017")

ggplot(counts_hgnc, aes(factor(0), counts)) + 
  facet_wrap(~hgnc_symbol, scales="free_y") +
  geom_violin() + 
  geom_boxplot(width=0.3) +
  geom_jitter(data=other, color="dark grey", size=1, alpha=0.9) +
  geom_point(data=highlight, color="red", size=2, alpha=0.9) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) + 
  labs(x="gene", 
       y="Expression (normalized counts)", 
       title="")
