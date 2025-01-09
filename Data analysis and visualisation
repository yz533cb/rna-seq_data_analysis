# Pre-processing RNA-seq data in R
library(tidyverse)

seqdata <- read_tsv("directory/counts/WT_gtr2_pib2.featureCounts", comment="#")

# Write a SampleInfo.txt file and combine with the seqdata to label the samples with ‘cell type’ and ‘status’
countdata <- seqdata %>%
  column_to_rownames("Geneid") %>% 
  # turn the geneid column into rownames
  rename_all(str_remove, ".sorted.bam") %>% 
  # remove the ".bam" from the column names
  select(sampleinfo$Sample) %>% 
  # keep sample columns using SampleInfo$Sample
  as.matrix()

# Filtering the genes
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]

# Library sizes bar plot
librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=20e6, lty=2)

# Count distribution boxplots
logcounts <- log2(countdata + 1)
statusCol <- match(sampleinfo$Status, c("diauxic_shift", "exponential_phase")) + 1
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)
abline(h=median(as.matrix(logcounts)), col="blue")

# Principal Component Analysis (PCA)
rlogcounts <- rlog(countdata)
pcDat <- prcomp(t(rlogcounts))
autoplot(pcDat,
         data = sampleinfo, 
         colour="CellType", 
         shape="Status",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample), box.padding = 0.8)

# Interactive MDS Plot with Glimma
library(Glimma)
glMDSPlot(rlogcounts, 
          labels = sampleinfo$Sample, 
          groups = sampleinfo[,c("CellType", "Status")], 
          folder = "directory/mds")

library(DESeq2)

# Create the design model formula
design <- as.formula(~ CellType + Status + CellType * Status)

# Set exponential_phase as the intercept
sampleinfo$Status <- factor(sampleinfo$Status, 
                            levels = c("exponential_phase", "diauxic_shift"))
sampleinfo$CellType <- factor(sampleinfo$CellType, 
                            levels = c("WT", "gtr2", "pib2"))
modelMatrix <- model.matrix(design, data = sampleinfo)

# Create the DESeqDataSet object
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design)
ddsObj <- DESeq(ddsObj.raw)

# Visualise the sample data distribution with PCA
vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup=c("Status", "CellType"))

# Visualise the comparison types in ddsObj
resultsNames(ddsObj)

# Generate a results table comparing gtr2∆ and WT at the diauxic shift phase (contrast given in the ddsObj)
resGtr2DS <- results(ddsObj, 
                  name="CellTypegtr2.Statusdiauxic_shift", 
                  alpha = 0.05)

# Generate a results table comparing gtr2∆ and pib2∆ (provide a contrast)
resGvP <- results(ddsObj,
                  contrast = c("CellType", "Gtr2", "Pib2")
                  alpha = 0.05)


library(apeglm)
library(ashr)

ddsShrinkGtr2DS <- lfcShrink(ddsObj, coef = "CellTypegtr2.Statusdiauxic_shift", type = "apeglm")

# for using the contract not provided in ddsObj
ddsShrinkGvP <- lfcShrink(dds = ddsObj, contrast = c("CellType", "gtr2", "pib2"), type = "ashr")

shrinkGtr2DS <- as.data.frame(ddsShrink) %>%
rownames_to_column("GeneID") %>%
filter(padj < 0.05 &
       !is.na(padj) & 
              abs(log2FoldChange) > 1 & 
       !is.na(GeneID))  %>%  
              rename(logFC=log2FoldChange, FDR=padj)

Yeast_annot_file <- read_tsv("directory/SCgff3_annotation.txt", comment= "#")
ShrinkGtr2DS_annot <- merge(Yeast_annot_file, shrinkGtr2DS, by="GeneID")

# p-value histogram
hist(shrinkGtr2DS$pvalue) 

# MA plot
plotMA(ddsShrinkGtr2DS, alpha=0.05)   

# ggplot2 with label on the pvalue top 10 genes
cutoff <- sort(shrinkGtr2DS$pvalue)[10]
shrinkGtr2DS <- shrinkGtr2DS %>% 
  mutate(TopGeneLabel = ifelse(pvalue <= cutoff, GeneID, ""))
ggplot(shrinkGtr2DS, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour = FDR < 0.05), shape = 20, size = 0.5) +
  geom_text(aes(label = TopGeneLabel)) +
  labs(x = "mean of normalised counts", y = "log fold change")

# Volcano plot
filtTab <- shrinkGtr2DS %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
ggplot(filtTab, aes(x = logFC, y = `-log10(FDR)`)) + 
  geom_point(aes(colour = FDR < 0.05), size = 1)

# Strip chart for gene expression
topgene <- filter(shrinkGtr2DS, GeneID == "YHR215W") # any gene of interest
geneID <- topgene$GeneID
plotCounts(ddsObj, gene = geneID, intgroup = c("CellType", "Status"),
           returnData = T) %>% 
  ggplot(aes(x = Status, y= log2(count))) +
  geom_point(aes(fill = Status), shape = 21, size = 2) +
  facet_wrap(~CellType) +
  expand_limits(y = 0)

# Interactive strip chart with Glimma
group <- str_remove_all(sampleinfo$Group, "[aeiou]")
de <- as.integer(shrinkGtr2DS$FDR <= 0.05 & !is.na(shrinkGtr2DS$FDR))
normCounts <- log2(counts(ddsObj))
glXYPlot(
  x = shrinkGtr2$logFC,
  y = -log10(shrinkGtr2$pvalue),
  xlab = "logFC",
  ylab = "FDR",
  main = "GTR2",
  counts = normCounts,
  groups = group,
  status = de,
  anno = shrinkGtr2[, c("GeneID", "Description")],
  folder = "directory/volcano"
)

library(ComplexHeatmap)
library(circlize)

sigGenes_Gtr2DS <- shrinkGtr2DS %>% 
  filter(FDR < 0.05 & !is.na(FDR) & 
         abs(logFC) > 1 & 
         !is.na(GeneID)) %>% 
  pull("GeneID")

plotDat_Gtr2DS <- vst(ddsObj)[sigGenes_Gtr2DS,] %>% 
  assay()
z.mat_Gtr2DS <- t(scale(t(plotDat_Gtr2DS), center=TRUE, scale=TRUE))

# Colour palette
myPalette <- c("green3", "black", "red3")
myRamp = colorRamp2(c(-2, 0, 2), myPalette)

# Annotate the heatmap
ha = HeatmapAnnotation(df = colData(ddsObj)[,c("CellType", "Status")], 
                       col = list("CellType" = c("WT" = "gold4", "gtr2" = "gold1", "pib2" = "yellow1"),
                                  "Status" = c("exponential_phase" = "deepskyblue1", "diauxic_shift" = "deepskyblue4")))

# Generate the heatmap
hmp_Gtr2DS <- Heatmap(z.mat_Gtr2DS,
                      name = "z-score",
                      col = myRamp,
                      cluster_rows = TRUE,   
                      row_dend_side = "right",
                      row_title_gp = gpar(fontsize = 5),
                      column_names_rot = 45,
                      cluster_columns = FALSE,
                      show_row_name = FALSE,
                      top_annotation = ha)
map <- draw(hmp_Gtr2DS, row_km = 12, row_km_repeats = 100)

roworder <- row_order(map)

# Combine multiple clusters
cluster1_grg <- append(roworder[["3"]], roworder[["4"]])
cluster2_grr <- append(roworder[["5"]], roworder[["9"]])
cluster3_rgg <- do.call(c, list(roworder[["6"]], roworder[["7"]], roworder[["10"]]))
cluster4_rgr <- roworder[["2"]]
cluster5_rrr <- append(roworder[["1"]], roworder[["8"]])

# Slice out clusters from shrink file
# Use cluster 1 as an example. The other clusters were analysed using the similar codes
cluster1grg_shrink <- shrinkGtr2DS %>%
  slice(roworder$`3`, roworder$`4`)

# Add annotation to the clusters
row_number_grg = length(cluster1_grg)
cluster1grg_annot = matrix(nrow = row_number_grg, ncol = 3, dimnames = list(NULL, c("GeneID", "GeneName", "Description")))

iteration = 0
for (i in cluster1_grg) {
     iteration = iteration + 1
     cluster1grg_annot[iteration, 1] <- tempGeneID <- rownames(z.mat_Gtr2DS)[i]
   
     iteration_1 = 0 
     while (iteration_1 < dim(sigGenes_Gtr2DS_annot)[1]) {
          iteration_1 = iteration_1 + 1
          if (tempGeneID == sigGenes_Gtr2DS_annot[iteration_1, 1]) {
            cluster1grg_annot[iteration, 2] <- sigGenes_Gtr2DS_annot[iteration_1, 2]
            cluster1grg_annot[iteration, 3] <- sigGenes_Gtr2DS_annot[iteration_1, 3]
          }
     }
}

# Export the gene names from one cluster (for online promoter analysis)
write(cluster1grg_annot[ ,1 ], “directory/GeneNames_cluster1grg.txt" )

library(goseq)

sigData_Gtr2DS <- as.integer( shrinkGtr2DS$FDR < 0.01 & !is.na(shrinkGtr2DS$FDR) & shrinkGtr2DS$logFC > 2)
names(sigData_Gtr2DS) <- shrinkGtr2DS$GeneID

# Fit the probability weighting function (PWF)
pwf_Gtr2DS <- nullp(sigData_Gtr2DS, "sacCer2", "ensGene", bias.data = shrinkGtr2DS$medianTxLength)

# Conduct GO enrichment analysis
goResults_Gtr2DS <- goseq(pwf_Gtr2DS, "sacCer2","ensGene", test.cats = c("GO:BP"), use_genes_without_cat = TRUE)
goResults_Gtr2DS %>% 
  # Plot the top 10 categories
  top_n( 10, wt = -over_represented_pvalue ) %>% 
  mutate( hitsPerc = numDEInCat*100/numInCat ) %>% 
  ggplot(aes(x = hitsPerc, 
             y = term, 
             colour = over_represented_pvalue, 
             size = numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs( x="Hits (%)", y="GO term", colour="p value", size="Count" )

# Receive the GO information for the GO accessions
# Use category 1 as an example
library(GO.db)
GOTERM[[goResults_Gtr2DS$category[1]]]

library(clusterProfiler)
library(pathview)

sigGenes_Gtr2DS <- shrinkGtr2DS %>% 
  filter(FDR < 0.05 & !is.na(FDR) & 
         abs(logFC) > 1 & 
         !is.na(GeneID)) %>% 
  pull(GeneID)
kk_Gtr2DS <- enrichKEGG(gene = sigGenes_Gtr2DS, organism = 'sce')
head(kk_Gtr2DS, n=10)

# Visualise a pathway
# Use sce03010 as an example
browseKEGG(kk_Gtr2DS_filter, 'sce03010')

# Pathways can be stored as a file
logFC_Gtr2DS <- sigGenes_Gtr2DS_annot$logFC
names(logFC_Gtr2DS) <- sigGenes_Gtr2DS_annot$GeneID
pathview(gene.data = logFC_Gtr2DS, 
         pathway.id = "sce00030", 
         gene.idtype = "KEGG",
         species = "sce", 
         kegg.native = T, same.layer=F,
         limit = list(gene=5, cpd=1))
