#Load libraries
#To install package use: BiocManager::install("packagename")
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.At.tair.db) #Arabidopsis thaliana genome annotation
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(pheatmap)
library(apeglm)
library(plotly)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)


######Import featureCounts output
# Import gene counts table
# - skip first row
# - make row names the gene identifiers
countdata = read.table("final_counts.txt", header = TRUE, as.is=T, skip = 1, row.names = 1)

# Remove Aligned.sortedByCoord.out.bam + '..' from column identifiers
colnames(countdata) = gsub("Aligned.sortedByCoord.out.bam", "", colnames(countdata), fixed = T)
colnames(countdata) = gsub("..", "", colnames(countdata), fixed = T)

# Remove "_R1_001" from column identifiers
colnames(countdata) = gsub("_R1_001", "", colnames(countdata), fixed = T)

# Remove first five columns
countdata = countdata[ ,c(-1:-5)]

# Make sure ID's are correct
head(countdata)


#######Import metadata csv file
# Import metadata file
metadata = read.csv("metadata.csv", as.is=T, header=T, row.names = 1)

# Reorder sampleID's to match featureCounts column order. 
metadata = metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
head(metadata)


######Make DESeq2 object from counts and metadata
# - heatmap : count dataframe (rounded to nearest integer)
# - colData : sample metadata in the dataframe (row names correspond to columns of heatmap)
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat_Col = DESeqDataSetFromMatrix(countData = round(countdata), colData = metadata, design = ~Group)
ddsMat_H4WT = DESeqDataSetFromMatrix(countData = round(countdata), colData = metadata, design = ~Group)

#Set the "reference" condition
ddsMat_Col$Group = relevel(ddsMat_Col$Group, ref = "Col")
ddsMat_H4WT$Group = relevel(ddsMat_H4WT$Group, ref = "H4WT")

#Calculate FPM
fpm_Col = fpm(ddsMat_Col)

# Find differential expressed genes
ddsMat_Col = DESeq(ddsMat_Col)
ddsMat_H4WT = DESeq(ddsMat_H4WT)

#Retrieve results table and make data frame
#Samples normalized to Col values
resultsNames(ddsMat_Col)
res_WTvsCol = results(ddsMat_Col, pAdjustMethod = "fdr", alpha = 0.05, name="Group_H4WT_vs_Col")
res_WTvsCol = as.data.frame(res_WTvsCol)
res_R17AvsCol = results(ddsMat_Col, pAdjustMethod = "fdr", alpha = 0.05, name="Group_R17A_vs_Col")
res_R17AvsCol = as.data.frame(res_R17AvsCol)
res_R35KvsCol = results(ddsMat_Col, pAdjustMethod = "fdr", alpha = 0.05, name="Group_R35K_vs_Col")
res_R35KvsCol = as.data.frame(res_R35KvsCol)
res_chrvsCol = results(ddsMat_Col, pAdjustMethod = "fdr", alpha = 0.05, name="Group_chr11chr17_vs_Col")
res_chrvsCol = as.data.frame(res_chrvsCol)
#Samples normalized to H4WT values
resultsNames(ddsMat_H4WT)
res_R17AvsH4WT = results(ddsMat_H4WT, pAdjustMethod = "fdr", alpha = 0.05, name="Group_R17A_vs_H4WT")
res_R17AvsH4WT = as.data.frame(res_R17AvsH4WT)
res_R35KvsH4WT = results(ddsMat_H4WT, pAdjustMethod = "fdr", alpha = 0.05, name="Group_R35K_vs_H4WT")
res_R35KvsH4WT = as.data.frame(res_R35KvsH4WT)
res_chrvsH4WT = results(ddsMat_H4WT, pAdjustMethod = "fdr", alpha = 0.05, name="Group_chr11chr17_vs_H4WT")
res_chrvsH4WT = as.data.frame(res_chrvsH4WT)

# Sort the results table by adjusted p-value and remove NA adjusted p-values
# Samples normalized to Col values
resordered_WTvsCol = data.frame(res_WTvsCol[order(res_WTvsCol$padj, na.last=NA),])
resordered_R17AvsCol = data.frame(res_R17AvsCol[order(res_R17AvsCol$padj, na.last=NA),])
resordered_R35KvsCol = data.frame(res_R35KvsCol[order(res_R35KvsCol$padj, na.last=NA),])
resordered_chrvsCol = data.frame(res_chrvsCol[order(res_chrvsCol$padj, na.last=NA),])
# Samples normalized to H4WT values
resordered_R17AvsH4WT = data.frame(res_R17AvsH4WT[order(res_R17AvsH4WT$padj, na.last=NA),])
resordered_R35KvsH4WT = data.frame(res_R35KvsH4WT[order(res_R35KvsH4WT$padj, na.last=NA),])
resordered_chrvsH4WT = data.frame(res_chrvsH4WT[order(res_chrvsH4WT$padj, na.last=NA),])

#How many adjusted p-values were less than 0.01?
sum(resordered_WTvsCol$padj < 0.01, na.rm=TRUE)
sum(resordered_R17AvsCol$padj < 0.01, na.rm=TRUE)
sum(resordered_R35KvsCol$padj < 0.01, na.rm=TRUE)
sum(resordered_chrvsCol$padj < 0.01, na.rm=TRUE)
sum(resordered_R17AvsH4WT$padj < 0.01, na.rm=TRUE)
sum(resordered_R35KvsH4WT$padj < 0.01, na.rm=TRUE)
sum(resordered_chrvsH4WT$padj < 0.01, na.rm=TRUE)

# Get genes with adjusted p-values less than 0.01 and log2FoldChange >= 1 and write to table
topDEgenes_WTvsCol <- rownames(resordered_WTvsCol[resordered_WTvsCol$padj<0.01 & resordered_WTvsCol$log2FoldChange >= 1,])
write.table(topDEgenes_WTvsCol, "topDEgenes_H4WTvsCol.txt", quote=F, col.names=F, row.names = F)
topDEgenes_R17AvsCol <- rownames(resordered_R17AvsCol[resordered_R17AvsCol$padj<0.01 & resordered_R17AvsCol$log2FoldChange >= 1,])
write.table(topDEgenes_R17AvsCol, "topDEgenes_R17AvsCol.txt", quote=F, col.names=F, row.names = F)
topDEgenes_R35KvsCol <- rownames(resordered_R35KvsCol[resordered_R35KvsCol$padj<0.01 & resordered_R35KvsCol$log2FoldChange >= 1,])
write.table(topDEgenes_R35KvsCol, "topDEgenes_R35KvsCol.txt", quote=F, col.names=F, row.names = F)
topDEgenes_chrvsCol <- rownames(resordered_chrvsCol[resordered_chrvsCol$padj<0.01 & resordered_chrvsCol$log2FoldChange >= 1,])
write.table(topDEgenes_chrvsCol, "topDEgenes_chrvsCol.txt", quote=F, col.names=F, row.names = F)
topDEgenes_R17AvsH4WT <- rownames(resordered_R17AvsH4WT[resordered_R17AvsH4WT$padj<0.01 & resordered_R17AvsH4WT$log2FoldChange >= 1,])
write.table(topDEgenes_R17AvsH4WT, "topDEgenes_R17AvsH4WT.txt", quote=F, col.names=F, row.names = F)
topDEgenes_R35KvsH4WT <- rownames(resordered_R35KvsH4WT[resordered_R35KvsH4WT$padj<0.01 & resordered_R35KvsH4WT$log2FoldChange >= 1,])
write.table(topDEgenes_R35KvsH4WT, "topDEgenes_R35KvsH4WT.txt", quote=F, col.names=F, row.names = F)
topDEgenes_chrvsH4WT <- rownames(resordered_chrvsH4WT[resordered_chrvsH4WT$padj<0.01 & resordered_chrvsH4WT$log2FoldChange >= 1,])
write.table(topDEgenes_chrvsH4WT, "topDEgenes_chrvsH4WT.txt", quote=F, col.names=F, row.names = F)




######Annotate gene symbols
#Function to annotate results and write to subsets of results to txt files
annotateAndPrint = function(results, ddsMat, output_name) {
  # Add gene full name
  results$description <- mapIds(x = org.At.tair.db, keys = row.names(results), column = "GENENAME", keytype = "TAIR", multiVals = "first")
  # Add gene symbol
  results$symbol <- row.names(results)
  # Add ENTREZ ID
  results$entrez <- mapIds(x = org.At.tair.db, keys = row.names(results), column = "ENTREZID", keytype = "TAIR", multiVals = "first")
  # Subset for only significant genes (p < 0.05)
  results_sig <- subset(results, padj < 0.05)
  # Write normalized gene counts to txt file
  write.table(x = as.data.frame(counts(ddsMat), normalized = T), file = paste(output_name, 'normalized_counts.txt'), sep = '\t', quote = F, col.names = NA)
  # Write significant normalized gene counts to a .txt file
  write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), file = paste(output_name, 'normalized_counts_significant.txt'), sep = '\t', quote = F, col.names = NA)
  # Write the annotated results table to a .txt file
  write.table(x = as.data.frame(results), file = paste(output_name, "results_gene_annotated.txt"), sep = '\t', quote = F, col.names = NA)
  # Write significant annotated results table to a .txt file
  write.table(x = as.data.frame(results_sig), file = paste(output_name, "results_gene_annotated_significant.txt"), sep = '\t', quote = F, col.names = NA)
  return(results_sig)
}   
results_sig_WTvsCol = annotateAndPrint(res_WTvsCol, ddsMat_Col, "H4WTvsCol")
results_sig_R17AvsCol = annotateAndPrint(res_R17AvsCol, ddsMat_Col, "R17AvsCol")
results_sig_R35KvsCol = annotateAndPrint(res_R35KvsCol, ddsMat_Col, "R35KvsCol")
results_sig_chrvsCol = annotateAndPrint(res_chrvsCol, ddsMat_Col, "chrvsCol")
results_sig_R17AvsH4WT = annotateAndPrint(res_R17AvsH4WT, ddsMat_H4WT, "R17AvsH4WT")
results_sig_R35KvsH4WT = annotateAndPrint(res_R35KvsH4WT, ddsMat_H4WT, "R35KvsH4WT")
results_sig_chrvsH4WT = annotateAndPrint(res_chrvsH4WT, ddsMat_H4WT, "chrvsH4WT")

####Determine up-regulated and down-regulated genes
R17AvsCol_up = subset(subset(res_R17AvsCol, padj<0.05),log2FoldChange>0)
R17AvsCol_down = subset(subset(res_R17AvsCol, padj<0.05),log2FoldChange<0)
chrvsCol_up = subset(subset(res_chrvsCol, padj<0.05),log2FoldChange>0)
chrvsCol_down = subset(subset(res_chrvsCol, padj<0.05),log2FoldChange<0)
R35KvsCol_up = subset(subset(res_R35KvsCol, padj<0.05),log2FoldChange>0)
R35KvsCol_down = subset(subset(res_R35KvsCol, padj<0.05),log2FoldChange<0)
H4WTvsCol_up = subset(subset(res_WTvsCol, padj<0.05),log2FoldChange>0)
H4WTvsCol_down = subset(subset(res_WTvsCol, padj<0.05),log2FoldChange<0)
R17AvsH4WT_up = subset(subset(res_R17AvsH4WT, padj<0.05),log2FoldChange>0)
R17AvsH4WT_down = subset(subset(res_R17AvsH4WT, padj<0.05),log2FoldChange<0)
chrvsH4WT_up = subset(subset(res_chrvsH4WT, padj<0.05),log2FoldChange>0)
chrvsH4WT_down = subset(subset(res_chrvsH4WT, padj<0.05),log2FoldChange<0)
R35KvsH4WT_up = subset(subset(res_R35KvsH4WT, padj<0.05),log2FoldChange>0)
R35KvsH4WT_down = subset(subset(res_R35KvsH4WT, padj<0.05),log2FoldChange<0)

R17AvsCol_upnames = rownames(R17AvsCol_up)
R17AvsCol_downnames = rownames(R17AvsCol_down)
chrvsCol_upnames = rownames(chrvsCol_up)
chrvsCol_downnames = rownames(chrvsCol_down)
R17Aandchr_up = intersect(R17AvsCol_upnames, chrvsCol_upnames)
R17Aandchr_down = intersect(R17AvsCol_downnames, chrvsCol_downnames)

H4WTvsCol_upnames = rownames(H4WTvsCol_up)
H4WTvsCol_downnames = rownames(H4WTvsCol_down)
R17AminusH4WT_upnames = setdiff(R17AvsCol_upnames, H4WTvsCol_upnames)
R17AminusH4WT_downnames = setdiff(R17AvsCol_downnames, H4WTvsCol_downnames)


#Plot Venn diagram of up-regulated and down-regulated genes
###Two-way
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
#function to make up-regulated Venn diagrams (2-way)
makeVenns = function(input1_up, input2_up, input1_down, input2_down, name1, name2) {
  venn.diagram(x=list(rownames(input1_up), rownames(input2_up)), category.names=c(name1, name2),
               filename=paste("Venn Diagrams/",name1,"and",name2,"Up-regulated.png"), output=T, lwd = 2, fill = c("lightblue","lightgreen"), 
               alpha=0.3, col=c("lightblue","lightgreen"), cat.pos=c(200,150), fontfamily = "sans", cat.fontfamily = "sans",
               cat.cex=2, cex=2, main.cex=2, sub.cex=2, main.fontfamily="sans", cat.dist = c(0.03, 0.05),
               sub=paste(name1,"and",name2), sub.fontfamily="sans", main="Up-regulated genes")
  venn.diagram(x=list(rownames(input1_down), rownames(input2_down)), category.names=c(name1, name2),
               filename=paste("Venn Diagrams/",name1,"and",name2,"Down-regulated.png"), output=T, lwd = 2, fill = c("pink","orange"),
               alpha=0.3, col=c("pink","orange"), cat.pos=c(200,150), fontfamily = "sans", cat.fontfamily = "sans",
               cat.cex=2, cex=2, main.cex=2, sub.cex=2, main.fontfamily="sans", cat.dist = c(0.03, 0.05),
               sub=paste(name1,"and",name2), sub.fontfamily="sans", main="Down-regulated genes")
}
#R17A and chr vs Col
makeVenns(R17AvsCol_up, chrvsCol_up, R17AvsCol_down, chrvsCol_down, "R17A", "chr11chr17")
#R35K and chr vs Col
makeVenns(R35KvsCol_up, chrvsCol_up, R35KvsCol_down, chrvsCol_down, "R35K", "chr11chr17")
#most significant R17A and chr vs Col
venn.diagram(x=list(rownames(R17AvsCol_mostsig), rownames(chrvsCol_mostsig)), category.names=c("R17A", "chr11chr17"),
             filename=paste("Venn Diagrams/Most significant R17A and chr DEG >8-fold.png"), output=T, lwd = 2, fill = c("mediumorchid2","slateblue1"), 
             alpha=0.3, col=c("mediumorchid2","slateblue1"), cat.pos=c(200,150), fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, sub.cex=2, main.fontfamily="sans", cat.dist = c(0.03, 0.05),
             sub="R17A and chr11chr17", sub.fontfamily="sans", main="Most Significant DEG")
#R17A (minus H4WT) and chr vs Col
venn.diagram(x=list(R17AminusH4WT_upnames, rownames(chrvsCol_up)), category.names=c("R17A", "chr11chr17"),
             filename=paste("Venn Diagrams/R17A (minus H4WT) and chr Up-regulated.png"), output=T, lwd = 2, fill = c("lightblue","lightgreen"), 
             alpha=0.3, col=c("lightblue","lightgreen"), cat.pos=c(200,150), fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, sub.cex=2, main.fontfamily="sans", cat.dist = c(0.03, 0.05),
             sub="R17A (minus H4WT) and chr11chr17", sub.fontfamily="sans", main="Up-regulated genes")
venn.diagram(x=list(R17AminusH4WT_downnames, rownames(chrvsCol_down)), category.names=c("R17A", "chr11chr17"),
             filename=paste("Venn Diagrams/R17A (minus H4WT) and chr Down-regulated.png"), output=T, lwd = 2, fill = c("pink","orange"), 
             alpha=0.3, col=c("pink","orange"), cat.pos=c(200,150), fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, sub.cex=2, main.fontfamily="sans", cat.dist = c(0.03, 0.05),
             sub="R17A (minus H4WT) and chr11chr17", sub.fontfamily="sans", main="Down-regulated genes")

#R17A, R35K, chr11chr17 and H4 WT compared to Col (up)
myCol <- brewer.pal(4, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(R17AvsCol_up), rownames(R35KvsCol_up), rownames(chrvsCol_up), rownames(H4WTvsCol_up)), 
             category.names=c("R17A","R35K","chr11chr17","H4 WT"),
             filename="Venn Diagrams/All Up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#R17A, R35K, chr11chr17 and H4 WT compared to Col (down)
venn.diagram(x=list(rownames(R17AvsCol_down), rownames(R35KvsCol_down), rownames(chrvsCol_down), rownames(H4WTvsCol_down)), 
             category.names=c("R17A","R35K","chr11chr17","H4 WT"),
             filename="Venn Diagrams/All Down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

#R17A, R35K, and chr11chr17 compared to Col (up)
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(R17AvsCol_up), rownames(R35KvsCol_up), rownames(chrvsCol_up)), 
             category.names=c("R17A","R35K","chr11chr17"),
             filename="Venn Diagrams/R17A-R35K-chr Up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#R17A, R35K, chr11chr17 and H4 WT compared to Col (down)
venn.diagram(x=list(rownames(R17AvsCol_down), rownames(R35KvsCol_down), rownames(chrvsCol_down)), 
             category.names=c("R17A","R35K","chr11chr17"),
             filename="Venn Diagrams/R17A-R35K-chr Down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

#R17A, H4WT, and chr11chr17 compared to Col (up)
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(R17AvsCol_up), rownames(chrvsCol_up), rownames(H4WTvsCol_up)), 
             category.names=c("R17A","chr11chr17","H4 WT"),
             filename="Venn Diagrams/R17A chr11 H4WT Up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#R17A, H4WT, chr11chr17 and H4 WT compared to Col (down)
venn.diagram(x=list(rownames(R17AvsCol_down), rownames(chrvsCol_down), rownames(H4WTvsCol_down)), 
             category.names=c("R17A","chr11chr17","H4 WT"),
             filename="Venn Diagrams/R17A chr11 H4WT Down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")



#####Plot Gene Expression Data

#PCA Plot
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat_Col, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 20) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 20 most variable genes") 
ggsave("PCA/top20.pdf", device=pdf())

#Heatmap
#all genes
mat <- assay(ddsMat_rlog[row.names(results_sig)])

###
res_sig_ordered_chrvsCol = subset(resordered_chrvsCol, padj <0.05)
mat <- assay(ddsMat_rlog[row.names(res_sig_ordered_chrvsCol)])[1:30,]
res_sig_ordered_R17AvsCol = subset(resordered_R17AvsCol, padj <0.05)
mat <- assay(ddsMat_rlog[row.names(res_sig_ordered_R17AvsCol)])[1:30,]
mat_chrvsCol <- assay(ddsMat_rlog[row.names(results_sig_chrvsCol)])
mat_R17AvsCol <- assay(ddsMat_rlog[row.names(results_sig_R17AvsCol)])
res_sig_ordered_R35KvsCol = subset(resordered_R35KvsCol, padj <0.05)
mat_R35KvsCol <- assay(ddsMat_rlog[row.names(results_sig_R35KvsCol)])
res_sig_ordered_WTvsCol = subset(resordered_WTvsCol, padj <0.05)
mat_WTvsCol <- assay(ddsMat_rlog[row.names(results_sig_WTvsCol)])


# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(Col = brewer.pal(5,"Set1")[1], chr11chr17 = brewer.pal(5,"Set1")[2], H4WT = brewer.pal(5,"Set1")[3], R17A = brewer.pal(5,"Set1")[4], R35K = brewer.pal(5,"Set1")[5]),
  Replicate = c(Rep1 = "gray", Rep2 = "black")
)

# Make Heatmap with pheatmap function
#Saves png to folder
plotheatmap = function(mat, name) {
  heatmap = pheatmap(mat = mat, 
           border_color = "white", 
           color= colorRampPalette(c("navy","white","firebrick3"))(30), 
           scale = "row", # Scale genes to Z-score (how many standard deviations)
           annotation_col = annotation_col, # Add multiple annotations to the samples
           annotation_colors = ann_colors,# Change the default colors of the annotations
           fontsize = 6.5, # Make fonts smaller
           cellwidth = 55, # Make the cells wider
           show_colnames = F, show_rownames = F, main = name)
}
heatmap_chrvsCol = plotheatmap(mat_chrvsCol, "DEG chr11 chr17")
heatmap_R17AvsCol = plotheatmap(mat_R17AvsCol, "DEG R17A")
heatmap_R35KvsCol = plotheatmap(mat_R35KvsCol, "DEG R35K")
heatmap_WTvsCol = plotheatmap(mat_WTvsCol, "DEG H4WT")

##Function to save heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  ggplot()
}

save_pheatmap_pdf(heatmap_chrvsCol, "heatmaps/ chr11 chr17 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_R17AvsCol, "heatmaps/ R17A DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_R35KvsCol, "heatmaps/ R35K DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_WTvsCol, "heatmaps/ H4 WT DEG heatmap.pdf", width=10, height=10)



##################
####Calculate correlation matrix using Pearson's correlation (graph Pearson's correlation plot using pheatmap)
#Normalized counts
normcounts = counts(ddsMat_Col, normalized=T)
labels = c()
for (i in 1:length(metadata$Group)) {
   labels[i] = paste(metadata$Group[i], metadata$Replicate[i])
}
colnames(normcounts) = labels
correlation_matrix = cor(normcounts, method="pearson")
corr_heatmap = pheatmap(mat = correlation_matrix, 
                        border_color = "black", 
                        color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                        show_colnames = T, show_rownames = T, main = "Correlation Matrix")
save_pheatmap_pdf(corr_heatmap, "Correlation/ Pearson's Correlation Matrix Normalized Counts.pdf", width=10, height=10)

#Log2fold changes
#By group
logmatrix = cbind(res_WTvsCol$log2FoldChange, res_R17AvsCol$log2FoldChange, res_R35KvsCol$log2FoldChange, res_chrvsCol$log2FoldChange)
colnames(logmatrix) = c("H4 WT", "R17A", "R35K", "chr11 chr17")
rownames(logmatrix) = rownames(res_WTvsCol)
correlation_log_matrix = cor(logmatrix, method="pearson", use="complete.obs")
corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                            border_color = "black", 
                            color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                            show_colnames = T, show_rownames = T, main = "Correlation Matrix")
save_pheatmap_pdf(corr_heatmap, "Correlation/ Pearson's Correlation Matrix Log2FoldChange.pdf", width=10, height=10)

#DEG genes in chr11chr17
logmatrix2 = logmatrix[res_chrvsCol$padj<0.05,]
colnames(logmatrix2) = c("H4WT-1", "H4WT-2", "R17A-1", "R17A-2", "R35K-1", "R35K-2", "chr11chr17-1", "chr11chr17-2")
correlation_log_matrix = cor(logmatrix2, method="pearson", use="complete.obs")
corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                        border_color = "black", 
                        color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                        show_colnames = T, show_rownames = T, main = "Correlation Matrix")
save_pheatmap_pdf(corr_heatmap, "Correlation/Pearson's Correlation Matrix Log2FoldChange chr11chr17 DEG.pdf", width=10, height=10)

#DEG genes in R17A
logmatrix2 = logmatrix[res_R17AvsCol$padj<0.05,]
colnames(logmatrix2) = c("H4WT-1", "H4WT-2", "R17A-1", "R17A-2", "R35K-1", "R35K-2", "chr11chr17-1", "chr11chr17-2")
correlation_log_matrix = cor(logmatrix2, method="pearson", use="complete.obs")
corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                        border_color = "black", 
                        color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                        show_colnames = T, show_rownames = T, main = "Correlation Matrix")
save_pheatmap_pdf(corr_heatmap, "Correlation/Pearson's Correlation Matrix Log2FoldChange R17A DEG.pdf", width=10, height=10)

#DEG genes in R35K
logmatrix2 = logmatrix[res_R35KvsCol$padj<0.05,]
colnames(logmatrix2) = c("H4WT-1", "H4WT-2", "R17A-1", "R17A-2", "R35K-1", "R35K-2", "chr11chr17-1", "chr11chr17-2")
correlation_log_matrix = cor(logmatrix2, method="pearson", use="complete.obs")
corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                        border_color = "black", 
                        color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                        show_colnames = T, show_rownames = T, main = "Correlation Matrix")
save_pheatmap_pdf(corr_heatmap, "Correlation/Pearson's Correlation Matrix Log2FoldChange R35K DEG.pdf", width=10, height=10)


#Volcano Plot
#Function to make volcano plot
volcanoplot = function(results, name) {
  # Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
  ## - Change pvalues to -log10 (1.3 = 0.05)
  data <- data.frame(gene = row.names(results), pval = -log10(results$padj), lfc = results$log2FoldChange)
  # Remove any rows that have NA as an entry
  data <- na.omit(data)
  # Color the points which are up or down
  ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
  ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
  data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))
  
  # Add ggplot2 layers
  vol +   
    ggtitle(label = paste("Volcano Plot",name), subtitle = "Colored by fold-change direction") +
    geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right") + # change the legend
    xlab(paste("log2(",name,")")) + # Change X-Axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  
  ggsave(paste("Volcano Plots/",name,".pdf"), device=pdf())
}
volcanoplot(res_WTvsCol, "H4WT v Col")
volcanoplot(res_R17AvsCol, "R17A v Col")
volcanoplot(res_R35KvsCol, "R35K v Col")
volcanoplot(res_chrvsCol, "chr11chr17 v Col")
volcanoplot(res_R17AvsH4WT, "R17A v H4WT")
volcanoplot(res_R35KvsH4WT, "R35K v H4WT")

#Function to make volcano plot with specific points labeled (takes in ordered results)
volcanoplot_labeled = function(results, name, points) {
  # Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
  ## - Change pvalues to -log10 (1.3 = 0.05)
  data <- data.frame(gene = row.names(results), pval = -log10(results$padj), lfc = results$log2FoldChange)
  # Remove any rows that have NA as an entry
  data <- na.omit(data)
  # Color the points which are up or down
  ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
  ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
  data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  # Add gene labels
  data <- data %>% mutate(genelabels = "")
  data$genelabels[1:points] <- rownames(results)[1:points]
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))
  
  # Add ggplot2 layers
  vol +   
    ggtitle(label = paste("Volcano Plot",name), subtitle = "Colored by fold-change direction") +
    geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
    geom_text_repel(aes(label = genelabels)) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right") + # change the legend
    xlab(paste("log2(",name,")")) + # Change X-Axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  
  ggsave(paste("Volcano Plots/",name,"labeled.pdf"), device=pdf())
}
volcanoplot_labeled(resordered_chrvsCol, "Chr vs Col", 10)
volcanoplot_labeled(resordered_R17AvsCol, "R17A vs Col", 10)
volcanoplot_labeled(resordered_R35KvsCol, "R35K vs Col", 10)
volcanoplot_labeled(resordered_WTvsCol, "H4WT vs Col", 10)



#####MA Plot
#MA Plot
make_MAplot_lfc = function(ddsMat, group) {
  resLFC <- lfcShrink(ddsMat, coef=group, type="apeglm")
  pdf(paste("MA Plots/",group,".pdf"))
  plotMA(resLFC, alpha=0.05, ylim = c(-5, 5), main=group)
  dev.off()
}
make_MAplot_lfc(ddsMat_Col, "Group_H4WT_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_R17A_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_R35K_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_chr11chr17_vs_Col")

make_MAplot = function(ddsMat, group) {
  res = results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05, name=group)
  pdf(paste("MA Plots/non-lfc",group,".pdf"))
  plotMA(res, alpha=0.05, ylim = c(-5, 5), main=group)
  dev.off()
}
make_MAplot(ddsMat_Col, "Group_H4WT_vs_Col")
make_MAplot(ddsMat_Col, "Group_R17A_vs_Col")
make_MAplot(ddsMat_Col, "Group_R35K_vs_Col")
make_MAplot(ddsMat_Col, "Group_chr11chr17_vs_Col")



########Single Gene Plots
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat_Col, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Single gene plots
plotCounts(dds = ddsMat_Col,gene = top_gene, intgroup = "Group", normalized = T, transform = T)

## Function to make single gene plots
singlegeneplot = function(genename, label, ddsMat) {
  d = plotCounts(dds = ddsMat, gene = genename, intgroup = "Group", normalized = T, transform = T, returnData = T)
  ggplot(d, aes(x = Group, y = count, color = Group)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    theme_bw() +
    ggtitle(label) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("Individual Gene Plots/",label,".pdf"), device=pdf())
}

####
#Flowering time genes
#Plot FT
singlegeneplot("AT1G65480", "FT", ddsMat_Col)
#Plot FLC
singlegeneplot("AT5G10140", "FLC", ddsMat_Col)
#Plot FLC4
singlegeneplot("AT5G65070","FLC4",ddsMat_Col)
#Plot FUL
singlegeneplot("AT5G60910", "FUL", ddsMat_Col)
#Plot SOC1
singlegeneplot("AT2G45660", "SOC1", ddsMat_Col)
#Plot SEP1
singlegeneplot("AT4G34190", "SEP1", ddsMat_Col)
#Plot SEP3
singlegeneplot("AT1G24260", "SEP3", ddsMat_Col)
#Plot CO
singlegeneplot("AT5G15840", "CO", ddsMat_Col)


#Multigene plot
multigeneplot = function(gene_list, label, normalized_counts, metadata, gene_names=NA) {
  gene_list <- normalized_counts %>% filter(gene %in% gene_list)
  # Add gene labels
  if (all(is.na(gene_names))==F) {
    gene_list <- gene_list %>% mutate(genenames = gene_names)
  }
  # Gathering the columns to have normalized counts to a single column
  gathered_list <- gene_list %>%
    gather(colnames(gene_list)[2:11], key = "sampleid", value = "normalized_counts")
  #combine metadata with normalized counts
  gathered_list <- inner_join(metadata, gathered_list)
  #Plot
  if (all(is.na(gene_names))==T) {
    ggplot(gathered_list) +
      geom_point(aes(x = gene, y = normalized_counts, color=Group)) +
      #scale_y_log10() +
      xlab("Genes") +
      ylab("Normalized Counts") +
      ggtitle(label) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  else {
    ggplot(gathered_list) +
      geom_point(aes(x = genenames, y = normalized_counts, color=Group)) +
      #scale_y_log10() +
      xlab("Genes") +
      ylab("Normalized Counts") +
      ggtitle(label) +
      theme_bw() +
      theme(text = element_text(size=20))
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  ggsave(paste("Individual Gene Plots/",label,".pdf"), device=pdf())
  dev.off()
}

#Make normalized counts data frame
normalized_counts = counts(ddsMat_Col, normalized=T)
normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
##Plot top 10 DEG in R17A
multigeneplot(rownames(resordered_R17AvsCol)[1:10], "Top 10 R17A DEG", normalized_counts, metadata)
##Plot top 20 DEG in R17A
multigeneplot(rownames(resordered_R17AvsCol)[1:20], "Top 20 R17A DEG", normalized_counts, metadata)
##Plot top 10 DEG in chr11chr17
multigeneplot(rownames(resordered_chrvsCol)[1:10], "Top 10 chr11chr17 DEG", normalized_counts, metadata)
##Plot top 20 DEG in chr11chr17
multigeneplot(rownames(resordered_chrvsCol)[1:20], "Top 20 chr11chr17 DEG", normalized_counts, metadata)
##Plot top 10 DEG in R35K
multigeneplot(rownames(resordered_R35KvsCol)[1:10], "Top 10 R35K DEG", normalized_counts, metadata)
##Plot top 20 DEG in R35K
multigeneplot(rownames(resordered_R35KvsCol)[1:20], "Top 20 R35K DEG", normalized_counts, metadata)
##Plot top 10 DEG in H4WT
multigeneplot(rownames(resordered_WTvsCol)[1:10], "Top 10 H4WT DEG", normalized_counts, metadata)
##Plot top 20 DEG in H4WT
multigeneplot(rownames(resordered_WTvsCol)[1:20], "Top 20 H4WT DEG", normalized_counts, metadata)
##Plot top 10 DEG shared between chr11 and R17A
multigeneplot(rownames(R17Aandchr_mostsig)[1:20], "Top 20 H4WT DEG", normalized_counts, metadata)

#######Flowering pathway analyses
#Make sublist of gene names and gene IDs
ftgenes = c("AT5G10140", "AT5G60910", "AT2G45660", "AT5G65070", "AT1G65480")
genenames_ft = c("FLC", "FUL", "SOC1", "FLC4", "FT")
ft_table = cbind(ftgenes, genenames_ft)
ft_table = ft_table[order(ft_table[,1]),]
#make gene plot
multigeneplot(ft_table[,1], "Top Flowering Genes", normalized_counts, metadata, gene_names=ft_table[,2])


##################
#Flowering time genes correlation
#read in flowering time gene data
florid_data = read.csv("/Users/emmacorcoran/Desktop/FLORID_all.csv")
#sort
florid_data = florid_data[order(florid_data$Gene.Details),]
#filter to only include flowering time genes
gene_list <- normalized_counts %>% filter(gene %in% florid_data$Gene.Details)
# Add gene labels
gene_list <- gene_list %>% mutate(genenames = florid_data$Short.Name)
# Gathering the columns to have normalized counts to a single column
gathered_list <- gene_list %>%
  gather(colnames(gene_list)[2:11], key = "sampleid", value = "normalized_counts")
#combine metadata with normalized counts
gathered_list <- inner_join(metadata, gathered_list)

#Make tibbles
#Make R17AvsCol data frame
R17AvsCol_counts <- res_R17AvsCol %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#Make chrvsCol data frame
R35KvsCol_counts <- res_R35KvsCol %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#Make chrvsCol data frame
chrvsCol_counts <- res_chrvsCol %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
#Make H4WTvsCol data frame
H4WTvsCol_counts <- res_WTvsCol %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#Filter rows to only include florid flowering genes and make matrix
R17A_flowering <- R17AvsCol_counts %>% filter(gene %in% florid_data$Gene.Details)
R17A_logs = R17A_flowering$log2FoldChange
R35K_flowering <- R35KvsCol_counts %>% filter(gene %in% florid_data$Gene.Details)
R35K_logs = R35K_flowering$log2FoldChange
chr_flowering <- chrvsCol_counts %>% filter(gene %in% florid_data$Gene.Details)
chr_logs = chr_flowering$log2FoldChange
H4WT_flowering <- H4WTvsCol_counts %>% filter(gene %in% florid_data$Gene.Details)
H4WT_logs = H4WT_flowering$log2FoldChange
flowering_matrix = cbind(R17A_flowering$log2FoldChange,R35K_flowering$log2FoldChange,
                         chr_flowering$log2FoldChange,H4WT_flowering$log2FoldChange)
rownames(flowering_matrix) = R17A_flowering$gene
colnames(flowering_matrix) = c("R17A","R35K","chr11chr17","H4WT")
flowering_matrix = data.frame(flowering_matrix)
flowering_matrix = na.omit(flowering_matrix)

#Make correlation scatterplot
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
#Initiate a ggplot
b <- ggplot(flowering_matrix, aes(x = R17A, y = chr11chr17, label=rownames(flowering_matrix))) + geom_point() +
  geom_smooth(method = lm, se = FALSE)
# Basic scatter plot + regression line
plotly_graph = ggplotly(b, tooltip = "all")
htmlwidgets::saveWidget(plotly_graph, "R17A_correlation.html")


#####All genes correlation
logchange_matrix = cbind(R17AvsCol_counts$log2FoldChange,R35KvsCol_counts$log2FoldChange,
                         chrvsCol_counts$log2FoldChange,H4WTvsCol_counts$log2FoldChange)
rownames(logchange_matrix) = H4WTvsCol_counts$gene
colnames(logchange_matrix) = c("R17A","R35K","chr11chr17","H4WT")
logchange_matrix = data.frame(logchange_matrix)
logchange_matrix = na.omit(logchange_matrix)

#Make correlation scatterplot
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
#Initiate a ggplot
b <- ggplot(logchange_matrix, aes(x = R17A, y = chr11chr17, label=rownames(logchange_matrix))) + geom_point() +
  geom_smooth(method = lm, se = FALSE)
# Basic scatter plot + regression line
plotly_graph = ggplotly(b, tooltip = "all")
htmlwidgets::saveWidget(plotly_graph, "R17A_all_correlation.html")

#####Regression plot
ggplotRegression <- function(fit){
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(lm(chr11chr17 ~ R17A, data = flowering_matrix))
ggsave("Regression/R17Avschr_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ R35K, data = flowering_matrix))
ggsave("Regression/R35Kvschr_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ H4WT, data = flowering_matrix))
ggsave("Regression/H4WTvschr_flowering_correlation.pdf")

ggplotRegression(lm(chr11chr17 ~ R17A, data = logchange_matrix))
ggsave("Regression/R17Avschr_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ R35K, data = logchange_matrix))
ggsave("Regression/R35Kvschr_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ H4WT, data = logchange_matrix))
ggsave("Regression/H4WTvschr_all_correlation.pdf")

###chr11chr17 DEG correlation
#Filter rows to only include chr11chr17 DEG and make matrix
R17A_chrDEG <- R17AvsCol_counts %>% filter(gene %in% rownames(results_sig_chrvsCol))
R35K_chrDEG <- R35KvsCol_counts %>% filter(gene %in% rownames(results_sig_chrvsCol))
chr_chrDEG <- chrvsCol_counts %>% filter(gene %in% rownames(results_sig_chrvsCol))
H4WT_chrDEG <- H4WTvsCol_counts %>% filter(gene %in% rownames(results_sig_chrvsCol))
chrDEG_matrix = cbind(R17A_chrDEG$log2FoldChange,R35K_chrDEG$log2FoldChange,
                         chr_chrDEG$log2FoldChange,H4WT_chrDEG$log2FoldChange)
rownames(chrDEG_matrix) = rownames(results_sig_chrvsCol)
colnames(chrDEG_matrix) = c("R17A","R35K","chr11chr17","H4WT")
chrDEG_matrix = data.frame(chrDEG_matrix)
chrDEG_matrix = na.omit(chrDEG_matrix)

#Initiate a ggplot
b <- ggplot(chrDEG_matrix, aes(x = R17A, y = chr11chr17, label=rownames(chrDEG_matrix))) + geom_point() +
  geom_smooth(method = lm, se = FALSE)
# Basic scatter plot + regression line
plotly_graph = ggplotly(b, tooltip = "all")
htmlwidgets::saveWidget(plotly_graph, "Regression/R17A_chr_correlation.html")
ggplotRegression(lm(chr11chr17 ~ R17A, data = chrDEG_matrix))
ggsave("Regression/R17Avschr_chrDEG_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ R35K, data = chrDEG_matrix))
ggsave("Regression/R35Kvschr_chrDEG_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ H4WT, data = chrDEG_matrix))
ggsave("Regression/H4WTvschr_chrDEG_correlation.pdf")


#####Finding Pathways for Differentially Expressed Genes
#http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
##Make list of overlapping entrez IDs
R17AvsCol_up_entrez = subset(results_sig_R17AvsCol, log2FoldChange>0)
chrvsCol_up_entrez = subset(results_sig_chrvsCol,log2FoldChange>0)
R17AvsCol_down_entrez = subset(results_sig_R17AvsCol, log2FoldChange<0)
chrvsCol_down_entrez = subset(results_sig_chrvsCol, log2FoldChange<0)
R17Aandchr_up_entrez = intersect(R17AvsCol_up_entrez$entrez, chrvsCol_up_entrez$entrez)
R17Aandchr_down_entrez = intersect(R17AvsCol_down_entrez$entrez, chrvsCol_down_entrez$entrez)


###Function to make barplots of KEGG Enrichment Pathways and GO Biological Pathways
pathway_plots = function (results, name) {
  # Create a matrix of gene log2 fold changes
  gene_matrix <- results$log2FoldChange
  # Add the entrezID's as names for each logFC entry
  names(gene_matrix) <- results$entrez
  #Enrich genes using KEGG database
  kegg_enrich <- enrichKEGG(gene = results$symbol, organism = 'ath', pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  # Plot results
  barplot(kegg_enrich, drop = TRUE, showCategory = 10, title = paste(name, "KEGG Enrichment Pathways"), font.size = 8)
  ggsave(paste("Pathways/", name, "KEGG Enrichment Pathways.pdf"), device=pdf())
  #Enrich genes using Gene Ontology
  go_enrich <- enrichGO(gene = names(gene_matrix), OrgDb = 'org.At.tair.db', readable = T,
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  # Plot results
  barplot(go_enrich,drop = TRUE, showCategory = 10, title = paste(name, "GO Biological Pathways"), font.size = 8)
  ggsave(paste("Pathways/",name,"GO Biological Pathways.pdf"), device=pdf())
}

#Generate barplots of pathways
pathway_plots(results_sig_R17AvsCol, "R17A")
pathway_plots(results_sig_R35KvsCol, "R35K")
pathway_plots(results_sig_chrvsCol, "chr11chr17")
pathway_plots(results_sig_H4WTvsCol, "H4 WT")
pathway_plots(results_sig_R17AvsH4WT, "R17A vs H4 WT")
pathway_plots(results_sig_R35KvsH4WT, "R35K vs H4 WT")
