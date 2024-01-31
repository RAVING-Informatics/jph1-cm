### 1) Read in outrider data set, add gene names, generate volcano plot

##Load packages
library(OUTRIDER)
library(annotables)

##Set working directory to whatever is suitable for your machine
setwd("/Users/00104561/Library/CloudStorage/OneDrive-TheUniversityofWesternAustralia/Research/Projects/RNAseq/DROP/Muscle_May2023/")

### Read in ods database results, and use bioMart to retrieve hgnc ids and  ###

##import ods.rds file and check the data structure
ods <- readRDS(file = "./processed_results/aberrant_expression/v38/outrider/MUSCLE/ods.Rds")

# remove versioning of gene IDs and merge with GRCh38 from annotables, obtaining HGNC symbols and other descriptions
geneIDs <- gsub("\\.[0-9]*(_[0-9]*)?.*$", "", rownames(ods))
map <- merge(data.table(ensgene=geneIDs), grch38, sort=FALSE,
             all.x=TRUE)[!duplicated(ensgene),] 

# set new gene names only if hgnc symbol is present
if(!"geneID" %in% colnames(mcols(ods))){
  mcols(ods)$ENSG <- geneIDs
  rownames(ods) <- map[,ifelse(
    is.na(symbol) | symbol == "" | duplicated(symbol), geneIDs, symbol)]
}

#checkpoint
rownames(ods)
colnames(assays(ods)$counts)

#save RDS
saveRDS(ods, file = "./ods-hgnc.Rds")

#Generate volcano plot

##Sample-based
plotVolcano(
  ods,
  basePlot=TRUE, 
  "D20-0017",
  col = c("gray", "firebrick"),
  main="JPH1",
  pch=1,
  padjCutoff = 0.05)

