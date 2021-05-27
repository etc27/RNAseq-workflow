rm(list=ls())

#https://github.com/twbattaglia/RNAseq-workflow

#Set working directory
setwd("/Users/emmacorcoran/Documents/R Scripts/RNA-seq/H4_ISWI_RNA-seq 3col")

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
library(apeglm)
library(plotly)
library(ggrepel)
library("VennDiagram")
library(tidyverse)


######Import featureCounts output
# Import gene counts table
# - skip first row
# - make row names the gene identifiers
countdata_0521 = read.table("/Users/emmacorcoran/Desktop/RNAseq_2021_05_21/4_final_counts/final_counts.txt", header = TRUE, as.is=T, skip = 1, row.names = 1)
countdata_1002 = read.table("/Users/emmacorcoran/Desktop/RNAseq_2020_10_02/5_final_counts/final_counts.txt", header = TRUE, as.is=T, skip = 1, row.names = 1)
# Remove first five columns
countdata_0521 = countdata_0521[ ,c(-1:-5)]
countdata_1002 = countdata_1002[ ,c(-1:-5)]
countdata = cbind(countdata_0521, countdata_1002)

#Reformat countdata
# Remove Aligned.sortedByCoord.out.bam + '..' from column identifiers
colnames(countdata) = gsub(".Aligned.sortedByCoord.out.bam", "", colnames(countdata), fixed = T)
colnames(countdata) = gsub("Aligned.sortedByCoord.out.bam", "", colnames(countdata), fixed = T)
colnames(countdata) = gsub("..", "", colnames(countdata), fixed = T)
# Remove "_R1_001" from column identifiers
colnames(countdata) = gsub("_R1_001", "", colnames(countdata), fixed = T)

# Make sure ID's are correct
head(countdata)

#remove chloroplast and mitochrondria reads - ATC... and ATM...
toremove = c()
for (i in 1:nrow(countdata)) {
  chrnum = substr(toString(rownames(countdata)[i]),3,3)
  if (chrnum == "C") {
    toremove = c(toremove,i)
  }
  else if(chrnum == "M") {
    toremove = c(toremove,i)
  }
}
fivechr_countdata = countdata[-toremove,]
countdata = fivechr_countdata
#Remove R35K, rlt12.2, and 68.18D.1
countdata = countdata[,-c(1,7,11,12)]
#Make sure countdata matrix looks correct
head(countdata)

#######Import metadata csv file
# Import metadata file
metadata = read.csv("metadata_3cols.csv", as.is=T, header=T, row.names = 1)

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

#Set the "reference" condition
ddsMat_Col$Group = relevel(ddsMat_Col$Group, ref = "Col")

#Calculate FPM
fpm_Col = fpm(ddsMat_Col)

# Find differential expressed genes
ddsMat_Col = DESeq(ddsMat_Col)

#Samples normalized to Col values
resultsNames(ddsMat_Col)

#Retrieve results table and make data frame
get_results = function(ddsMat, group_name){
  output_res = results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05, name=group_name)
  output_res = as.data.frame(output_res)
  resordered = data.frame(output_res[order(output_res$padj, na.last=NA),])
  print(paste("Number of significant DEGs:", sum(resordered$padj < 0.05 & abs(resordered$log2FoldChange) > 1, na.rm=TRUE)))
  topDEgenes <- rownames(resordered[resordered$padj<0.05 & abs(resordered$log2FoldChange)>1,])
  write.table(topDEgenes, paste0("topDEgenes/",group_name), quote=F, col.names=F, row.names = F)
  return(output_res)
}
res_68.18DvsCol = get_results(ddsMat_Col, "Group_68.18D_vs_Col")
res_68.21CvsCol = get_results(ddsMat_Col, "Group_68.21C_vs_Col")
res_13.6vsCol = get_results(ddsMat_Col, "Group_13.6_vs_Col")
res_13.9vsCol = get_results(ddsMat_Col, "Group_13.9_vs_Col")
res_arid5vsCol = get_results(ddsMat_Col, "Group_arid5_vs_Col")
res_rlt12vsCol = get_results(ddsMat_Col, "Group_rlt12_vs_Col")
res_chrvsCol = get_results(ddsMat_Col, "Group_chr11chr17_vs_Col")

#Annotate gene symbols
#Function to annotate results and write to subsets of results to txt files
annotateAndPrint = function(results, ddsMat, output_name) {
  # Add gene full name
  results$description <- mapIds(x = org.At.tair.db, keys = row.names(results), column = "GENENAME", keytype = "TAIR", multiVals = "first")
  # Add gene symbol
  results$symbol <- row.names(results)
  # Add ENTREZ ID
  results$entrez <- mapIds(x = org.At.tair.db, keys = row.names(results), column = "ENTREZID", keytype = "TAIR", multiVals = "first")
  # Subset for only significant genes (p < 0.05)
  results_sig <- subset(results, padj < 0.05 & abs(results$log2FoldChange) > 1)
  # Write normalized gene counts to txt file
  write.table(x = as.data.frame(counts(ddsMat), normalized = T), file = paste0("normalized_counts/", output_name, '_normalized_counts.txt'), sep = '\t', quote = F, col.names = NA)
  # Write significant normalized gene counts to a .txt file
  write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), file = paste0("normalized_counts/",output_name, '_normalized_counts_significant.txt'), sep = '\t', quote = F, col.names = NA)
  # Write the annotated results table to a .txt file
  write.table(x = as.data.frame(results), file = paste0("annotations/", output_name, "_results_gene_annotated.txt"), sep = '\t', quote = F, col.names = NA)
  # Write significant annotated results table to a .txt file
  write.table(x = as.data.frame(results_sig), file = paste0("annotations/", output_name, "_results_gene_annotated_significant.txt"), sep = '\t', quote = F, col.names = NA)
  return(results_sig)
}
results_sig_68.18DvsCol = annotateAndPrint(res_68.18DvsCol, ddsMat_Col, "68.18DvsCol")
results_sig_68.21CvsCol = annotateAndPrint(res_68.21CvsCol, ddsMat_Col, "68.21CvsCol")
results_sig_13.6vsCol = annotateAndPrint(res_13.6vsCol, ddsMat_Col, "13.6vsCol")
results_sig_13.9vsCol = annotateAndPrint(res_13.9vsCol, ddsMat_Col, "13.9vsCol")
results_sig_arid5vsCol = annotateAndPrint(res_arid5vsCol, ddsMat_Col, "arid5vsCol")
results_sig_rlt12vsCol = annotateAndPrint(res_rlt12vsCol, ddsMat_Col, "rlt12vsCol")
results_sig_chrvsCol = annotateAndPrint(res_chrvsCol, ddsMat_Col, "chrvsCol")

#Determine up-regulated and down-regulated genes
up_genes = function(results, sample){
  up_genes = subset(results, padj<0.05&log2FoldChange>1)
  upnames = rownames(up_genes)
  write.table(x = upnames, file = paste0("Up_and_Down/",sample,"_upreg.txt"), sep = '\t', quote = F, col.names = F, row.names=F)
  print(paste("Number of up-regulated genes:", nrow(up_genes)))
  return(up_genes)
}
down_genes = function(results, sample){
  down_genes = subset(results, padj<0.05&log2FoldChange<(-1))
  downnames = rownames(down_genes)
  write.table(x = downnames, file = paste0("Up_and_Down/",sample,"_downreg.txt"), sep = '\t', quote = F, col.names = F, row.names=F)
  print(paste("Number of down-regulated genes:", nrow(down_genes)))
  return(down_genes)
}
s68.18DvsCol_up = up_genes(res_68.18DvsCol, "68.18D")
s68.18DvsCol_down = down_genes(res_68.18DvsCol, "68.18D")
s68.21CvsCol_up = up_genes(res_68.21CvsCol, "68.21C")
s68.21CvsCol_down = down_genes(res_68.21CvsCol, "68.21C")
s13.6vsCol_up = up_genes(res_13.6vsCol, "13.6")
s13.6vsCol_down = down_genes(res_13.6vsCol, "13.6")
s13.9vsCol_up = up_genes(res_13.9vsCol, "13.9")
s13.9vsCol_down = down_genes(res_13.9vsCol, "13.9")
rlt12vsCol_up = up_genes(res_rlt12vsCol, "rlt12")
rlt12vsCol_down = down_genes(res_rlt12vsCol, "rlt12")
arid5vsCol_up = up_genes(res_arid5vsCol, "arid5")
arid5vsCol_down = down_genes(res_arid5vsCol, "arid5")
chrvsCol_up = up_genes(res_chrvsCol, "chr11chr17")
chrvsCol_down = down_genes(res_chrvsCol, "chr11chr17")

###hypergeometric test for DEG overlap
#all DEG (13.6 and arid5 overlap)
arid5and13.6_up = intersect(rownames(arid5vsCol_up), rownames(s13.6vsCol_up))
arid5and13.6_down = intersect(rownames(arid5vsCol_down), rownames(s13.6vsCol_down))
total_genes = nrow(ddsMat_Col)
overlap = length(c(arid5and13.6_up, arid5and13.6_down))
group1 = nrow(results_sig_arid5vsCol)
group2 = nrow(results_sig_13.6vsCol)
1-phyper(overlap, group1, total_genes-group1, group2)

#up-regulated
total_genes = nrow(ddsMat_Col)
overlap = length(arid5and13.6_up)
group1 = nrow(arid5vsCol_up)
group2 = nrow(s13.6vsCol_up)
1-phyper(overlap, group1, total_genes-group1, group2)

#down-regulated
total_genes = nrow(ddsMat_Col)
overlap = length(arid5and13.6_down)
group1 = nrow(arid5vsCol_down)
group2 = nrow(s13.6vsCol_down)
1-phyper(overlap, group1, total_genes-group1, group2)

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
#13.6 and arid5
makeVenns(s13.6vsCol_up, arid5vsCol_up, s13.6vsCol_down, arid5vsCol_down, "13.6", "arid5")
#13.6 and rlt12
makeVenns(s13.6vsCol_up, rlt12vsCol_up, s13.6vsCol_down, rlt12vsCol_down, "13.6", "rlt12")
#13.6 and chr11chr17
makeVenns(s13.6vsCol_up, chrvsCol_up, s13.6vsCol_down, chrvsCol_down, "13.6", "chr11chr17")
#arid5 and rlt12
makeVenns(arid5vsCol_up, rlt12vsCol_up, arid5vsCol_down, rlt12vsCol_down, "arid5", "rlt12")
#13.6 and 68.18D
makeVenns(s13.6vsCol_up, s68.18DvsCol_up, s13.6vsCol_down, s68.18DvsCol_down, "13.6", "68.18D")
#rlt12 and 68.18D
makeVenns(rlt12vsCol_up, s68.18DvsCol_up, rlt12vsCol_down, s68.18DvsCol_down, "rlt12", "68.18D")
#arid5 and 68.18D
makeVenns(arid5vsCol_up, s68.18DvsCol_up, arid5vsCol_down, s68.18DvsCol_down, "arid5", "68.18D")
#13.9 and arid5
makeVenns(s13.9vsCol_up, arid5vsCol_up, s13.9vsCol_down, arid5vsCol_down, "13.9", "arid5")
#13.9 and rlt12
makeVenns(s13.9vsCol_up, rlt12vsCol_up, s13.9vsCol_down, rlt12vsCol_down, "13.9", "rlt12")
#13.9 and rlt12
makeVenns(s13.9vsCol_up, s13.6vsCol_up, s13.9vsCol_down, s13.6vsCol_down, "13.9", "13.6")

###Venn diagrams with 5 samples - difficult to visualize
#all samples compared to Col (up)
myCol <- brewer.pal(5, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(s13.6vsCol_up), rownames(s13.9vsCol_up), rownames(arid5vsCol_up), rownames(rlt12vsCol_up), rownames(chrvsCol_up)), 
             category.names=c("13.6","13.9", "arid5", "rlt12", "chr11chr17"),
             filename="Venn Diagrams/All Up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#all samples compared to Col (down)
venn.diagram(x=list(rownames(s13.6vsCol_down), rownames(s13.9vsCol_down), rownames(arid5vsCol_down), rownames(rlt12vsCol_down), rownames(chrvsCol_down)), 
             category.names=c("13.6","13.9", "arid5", "rlt12", "chr11chr17"),
             filename="Venn Diagrams/All Down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

###Venn diagrams with 3 samples
#ISWI subunit mutants
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(arid5vsCol_up), rownames(rlt12vsCol_up), rownames(chrvsCol_up)), 
             category.names=c("arid5", "rlt12", "chr11chr17"),
             filename="Venn Diagrams/ISWI_up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#all samples compared to Col (down)
venn.diagram(x=list(rownames(arid5vsCol_down), rownames(rlt12vsCol_down), rownames(chrvsCol_down)), 
             category.names=c("arid5", "rlt12", "chr11chr17"),
             filename="Venn Diagrams/ISWI_down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

#chr11chr17 with 13.6 and 13.9
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(s13.6vsCol_up), rownames(s13.9vsCol_up), rownames(chrvsCol_up)), 
             category.names=c("13.6", "13.9", "chr11chr17"),
             filename="Venn Diagrams/13_chr_up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#all samples compared to Col (down)
venn.diagram(x=list(rownames(s13.6vsCol_down), rownames(s13.9vsCol_down), rownames(chrvsCol_down)), 
             category.names=c("13.6", "13.9", "chr11chr17"),
             filename="Venn Diagrams/13_chr_down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

#arid5 with 13.6 and 13.9
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(s13.6vsCol_up), rownames(s13.9vsCol_up), rownames(arid5vsCol_up)), 
             category.names=c("13.6", "13.9", "arid5"),
             filename="Venn Diagrams/13_arid5_up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#all samples compared to Col (down)
venn.diagram(x=list(rownames(s13.6vsCol_down), rownames(s13.9vsCol_down), rownames(arid5vsCol_down)), 
             category.names=c("13.6", "13.9", "arid5"),
             filename="Venn Diagrams/13_arid5_down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

#rlt12 with 13.6 and 13.9
myCol <- brewer.pal(3, "Pastel2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") #suppress log files saving
venn.diagram(x=list(rownames(s13.6vsCol_up), rownames(s13.9vsCol_up), rownames(rlt12vsCol_up)), 
             category.names=c("13.6", "13.9", "rlt12"),
             filename="Venn Diagrams/13_rlt12_up.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Up-regulated genes")
#all samples compared to Col (down)
venn.diagram(x=list(rownames(s13.6vsCol_down), rownames(s13.9vsCol_down), rownames(rlt12vsCol_down)), 
             category.names=c("13.6", "13.9", "rlt12"),
             filename="Venn Diagrams/13_rlt12_down.png", output=T, lwd = 2, fill = myCol,
             alpha=0.3, col=myCol, fontfamily = "sans", cat.fontfamily = "sans",
             cat.cex=2, cex=2, main.cex=2, main.fontfamily="sans",
             main="Down-regulated genes")

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


#Significant genes
sig_genes = function(results){
  output_res = subset(results, padj<0.05&abs(log2FoldChange)>1)
  resordered = data.frame(output_res[order(output_res$padj, na.last=NA),])
  print(paste("Number of significant genes:", nrow(resordered)))
  return(resordered)
}
#Make matrix
make_mat = function(results, rlog){
  mat = assay(rlog[row.names(results)])
  return(mat)
}
res_sig_ordered_68.18DvsCol = sig_genes(res_68.18DvsCol)
mat_68.18DvsCol = make_mat(res_sig_ordered_68.18DvsCol, ddsMat_rlog)
res_sig_ordered_68.21CvsCol = sig_genes(res_68.21CvsCol)
mat_68.21CvsCol = make_mat(res_sig_ordered_68.21CvsCol, ddsMat_rlog)
res_sig_ordered_13.6vsCol = sig_genes(res_13.6vsCol)
mat_13.6vsCol = make_mat(res_sig_ordered_13.6vsCol, ddsMat_rlog)
res_sig_ordered_13.9vsCol = sig_genes(res_13.9vsCol)
mat_13.9vsCol = make_mat(res_sig_ordered_13.9vsCol, ddsMat_rlog)
res_sig_ordered_arid5vsCol = sig_genes(res_arid5vsCol)
mat_arid5vsCol = make_mat(res_sig_ordered_arid5vsCol, ddsMat_rlog)
res_sig_ordered_rlt12vsCol = sig_genes(res_rlt12vsCol)
mat_rlt12vsCol = make_mat(res_sig_ordered_rlt12vsCol, ddsMat_rlog)
res_sig_ordered_chrvsCol = sig_genes(res_chrvsCol)
mat_chrvsCol = make_mat(res_sig_ordered_chrvsCol, ddsMat_rlog)

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(Col = brewer.pal(8,"Set1")[1], rlt12 = brewer.pal(8,"Set1")[2], arid5 = brewer.pal(8,"Set1")[3], "13.6" = brewer.pal(8,"Set1")[4], "13.9" = brewer.pal(8,"Set1")[5], "68.18D" = brewer.pal(8,"Set1")[6], "68.21C" = brewer.pal(8,"Set1")[7], "chr11chr17" = brewer.pal(8,"Set1")[8]),
  Replicate = c(Rep1 = "gray", Rep2 = "black")
)

# Make Heatmap with pheatmap function.
# Individual replicates
plotheatmap = function(mat, name) {
  heatmap = pheatmap(mat = mat, 
           border_color = "white", 
           color= colorRampPalette(c("navy","white","firebrick3"))(30), 
           scale = "row", # Scale genes to Z-score (how many standard deviations)
           annotation_col = metadata[,1:2], # Add multiple annotations to the samples
           #annotation_col = annotation_col, # Add multiple annotations to the samples
           annotation_colors = ann_colors,# Change the default colors of the annotations
           fontsize = 6.5, # Make fonts smaller
           cellwidth = 40, # Make the cells wider
           show_colnames = F, show_rownames = F, main = name)
}
heatmap_13.6vsCol = plotheatmap(mat_13.6vsCol, "DEG 13.6")
heatmap_13.9vsCol = plotheatmap(mat_13.9vsCol, "DEG 13.9")
heatmap_rlt12vsCol = plotheatmap(mat_rlt12vsCol, "DEG rlt12")
heatmap_arid5vsCol = plotheatmap(mat_arid5vsCol, "DEG arid5")
heatmap_chrvsCol = plotheatmap(mat_chrvsCol, "DEG chr")

##Function to save heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap_13.6vsCol, "heatmaps/13.6 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_13.9vsCol, "heatmaps/13.9 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_rlt12vsCol, "heatmaps/rlt12 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_arid5vsCol, "heatmaps/arid5 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(heatmap_chrvsCol, "heatmaps/chr DEG heatmap.pdf", width=10, height=10)

# Merged replicates
make_merged_mat = function(mat_to_merge){
  merged_mat = cbind(rowMeans(mat_to_merge[,c("rlt12.1","rlt12.1")]), rowMeans(mat_to_merge[,c("Col.1","Col.2","Col.3")]), 
                       rowMeans(mat_to_merge[,c("arid5.1","arid5.2")]), rowMeans(mat_to_merge[,c("chr11chr17.1","chr11chr17.2")]),
                       rowMeans(mat_to_merge[,c("13.6.1","13.6.2")]), rowMeans(mat_to_merge[,c("13.9.1","13.9.2")]),
                       rowMeans(mat_to_merge[,c("68.21C.1","68.21C.2")]), rowMeans(mat_to_merge[,c("68.18D.2","68.18D.2")]))
  colnames(merged_mat) = c("rlt12","Col","arid5","chr11chr17","13.6","13.9","68.21C","68.18D")
  return(merged_mat)
}
merged_mat_13.6vsCol = make_merged_mat(mat_13.6vsCol)
merged_mat_13.9vsCol = make_merged_mat(mat_13.9vsCol)
merged_mat_68.18DvsCol = make_merged_mat(mat_68.18DvsCol)
merged_mat_68.21CvsCol = make_merged_mat(mat_68.21CvsCol)
merged_mat_arid5vsCol = make_merged_mat(mat_arid5vsCol)
merged_mat_rlt12vsCol = make_merged_mat(mat_rlt12vsCol)
merged_mat_chrvsCol = make_merged_mat(mat_chrvsCol)
#Make heatmaps
plot_merged_heatmap = function(mat, name) {
  heatmap = pheatmap(mat = mat, 
                     border_color = "white", 
                     color= colorRampPalette(c("navy","white","firebrick3"))(30), 
                     scale = "row", # Scale genes to Z-score (how many standard deviations)
                     #annotation_col = ann_colors$Group, # Add multiple annotations to the samples
                     annotation_colors = ann_colors,# Change the default colors of the annotations
                     fontsize = 6.5, # Make fonts smaller
                     cellwidth = 55, # Make the cells wider
                     show_colnames = T, show_rownames = F, main = name)
}
merged_heatmap_13.6vsCol = plot_merged_heatmap(merged_mat_13.6vsCol, "DEG 13.6")
merged_heatmap_13.9vsCol = plot_merged_heatmap(merged_mat_13.9vsCol, "DEG 13.9")
merged_heatmap_rlt12vsCol = plot_merged_heatmap(merged_mat_rlt12vsCol, "DEG rlt12")
merged_heatmap_arid5vsCol = plot_merged_heatmap(merged_mat_arid5vsCol, "DEG arid5")
merged_heatmap_chrvsCol = plot_merged_heatmap(merged_mat_chrvsCol, "DEG chr")
#Save heatmaps
save_pheatmap_pdf(merged_heatmap_13.6vsCol, "merged heatmaps/13.6 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(merged_heatmap_13.9vsCol, "merged heatmaps/13.9 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(merged_heatmap_rlt12vsCol, "merged heatmaps/rlt12 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(merged_heatmap_arid5vsCol, "merged heatmaps/arid5 DEG heatmap.pdf", width=10, height=10)
save_pheatmap_pdf(merged_heatmap_chrvsCol, "merged heatmaps/chr DEG heatmap.pdf", width=10, height=10)

#Calculate correlation matrix using Pearson's correlation
#Log2fold changes
#By group
logmatrix = cbind(res_68.18DvsCol$log2FoldChange, res_68.21CvsCol$log2FoldChange, res_13.6vsCol$log2FoldChange, res_13.9vsCol$log2FoldChange, res_arid5vsCol$log2FoldChange, res_rlt12vsCol$log2FoldChange, res_chrvsCol$log2FoldChange)
colnames(logmatrix) = c("68.18D", "68.21C", "13.6", "13.9", "arid5", "rlt12", "chr11chr17")
rownames(logmatrix) = rownames(res_68.18DvsCol)

corr_matrix = function(logmatrix, name){
  correlation_log_matrix = cor(logmatrix, method="pearson", use="complete.obs")
  corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                          border_color = "black", 
                          color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                          show_colnames = T, show_rownames = T, main = "Correlation Matrix")
  save_pheatmap_pdf(corr_heatmap, paste0("Correlation/Pearson/",name,".pdf"), width=10, height=10)
  correlation_log_matrix = cor(logmatrix, method="spearman", use="complete.obs")
  corr_heatmap = pheatmap(mat = correlation_log_matrix, 
                          border_color = "black", 
                          color= colorRampPalette(c("white","navy"))(30), cluster_rows=T, cluster_cols=T,
                          show_colnames = T, show_rownames = T, main = "Correlation Matrix")
  save_pheatmap_pdf(corr_heatmap, paste0("Correlation/Spearman/",name,".pdf"), width=10, height=10)
}
#all genes
corr_matrix(logmatrix, "all_log2fold")
#13.6 DEGs
logmatrix_13.6 = na.omit(logmatrix[res_13.6vsCol$padj<0.05 & abs(res_13.6vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_13.6, "13.6")
#13.9 DEGs
logmatrix_13.9 = na.omit(logmatrix[res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_13.9, "13.9")
#13.6 and 13.9 DEGs
logmatrix_13.9 = na.omit(logmatrix[res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_13.9, "13.9")
#13.6 and 13.9
logmatrix_13 = na.omit(logmatrix[res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1 | res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_13, "13.6and13.9")
#chr11chr17
logmatrix_chr = na.omit(logmatrix[res_chrvsCol$padj<0.05 & abs(res_chrvsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_chr, "chr11chr17")
#rlt12
logmatrix_rlt12 = na.omit(logmatrix[res_rlt12vsCol$padj<0.05 & abs(res_rlt12vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_rlt12, "rlt12")
#arid5
logmatrix_arid5 = na.omit(logmatrix[res_arid5vsCol$padj<0.05 & abs(res_arid5vsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_arid5, "arid5")
#All DEGS
logmatrix_all = na.omit(logmatrix[res_68.21CvsCol$padj<0.05 & abs(res_68.21CvsCol$log2FoldChange)>1 | res_68.18DvsCol$padj<0.05 & abs(res_68.18DvsCol$log2FoldChange)>1 |
                            res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1 | res_13.6vsCol$padj<0.05 & abs(res_13.6vsCol$log2FoldChange)>1 |
                            res_arid5vsCol$padj<0.05 & abs(res_arid5vsCol$log2FoldChange)>1 | res_rlt12vsCol$padj<0.05 & abs(res_rlt12vsCol$log2FoldChange)>1 |
                            res_chrvsCol$padj<0.05 & abs(res_chrvsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_all, "all")
#All DEGS minus 68
logmatrix_allminus68 = na.omit(logmatrix[res_13.9vsCol$padj<0.05 & abs(res_13.9vsCol$log2FoldChange)>1 | res_13.6vsCol$padj<0.05 & abs(res_13.6vsCol$log2FoldChange)>1 |
                            res_arid5vsCol$padj<0.05 & abs(res_arid5vsCol$log2FoldChange)>1 | res_rlt12vsCol$padj<0.05 & abs(res_rlt12vsCol$log2FoldChange)>1 |
                            res_chrvsCol$padj<0.05 & abs(res_chrvsCol$log2FoldChange)>1,])
corr_matrix(logmatrix_allminus68, "allminus68")

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
  
  ggsave(paste0("Volcano Plots/",name,".pdf"), device=pdf())
  dev.off()
}
volcanoplot(res_68.18DvsCol, "68.18D v Col")
volcanoplot(res_68.21CvsCol, "68.21C v Col")
volcanoplot(res_13.6vsCol, "13.6 v Col")
volcanoplot(res_13.9vsCol, "13.9 v Col")
volcanoplot(res_arid5vsCol, "arid5 v Col")
volcanoplot(res_rlt12vsCol, "rlt12 v Col")
volcanoplot(res_chrvsCol, "chr11chr17 v Col")

#MA Plot
make_MAplot_lfc = function(ddsMat, group) {
  resLFC <- lfcShrink(ddsMat, coef=group, type="apeglm")
  pdf(paste0("MA Plots/",group,".pdf"))
  plotMA(resLFC, alpha=0.05, ylim = c(-5, 5), main=group)
  dev.off()
}
make_MAplot_lfc(ddsMat_Col, "Group_68.18D_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_68.21C_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_13.6_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_13.9_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_arid5_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_rlt12_vs_Col")
make_MAplot_lfc(ddsMat_Col, "Group_chr11chr17_vs_Col")

#Single Gene Plots
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat_Col, blind = FALSE)

## Function to make single gene plots
singlegeneplot = function(genename, label, ddsMat) {
  d = plotCounts(dds = ddsMat, gene = genename, intgroup = "Group", normalized = T, transform = T, returnData = T)
  ggplot(d, aes(x = Group, y = count, color = Group)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    theme_bw() +
    ggtitle(label) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Individual Gene Plots/",label,".pdf"), device=pdf())
}

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

#Plot H4 (At3g53730)
singlegeneplot("AT3G53730", "AT3G53730 (H4)", ddsMat_Col)
#Plot H4 (At2g28740)
singlegeneplot("AT2G28740", "AT2G28740 (H4)", ddsMat_Col)
#Plot chr11
singlegeneplot("AT3G06400", "CHR11", ddsMat_Col)
#Plot chr17
singlegeneplot("AT5G18620", "CHR17", ddsMat_Col)
#Plot REPRESSOR OF SILENCING1
singlegeneplot("AT2G36490", "ROS1", ddsMat_Col)

#DNA damage analyses
singlegeneplot("AT4G21070", "BRCA1", ddsMat_Col)
singlegeneplot("AT2G31970", "RAD50", ddsMat_Col)
singlegeneplot("AT1G20750", "UVH6", ddsMat_Col)
singlegeneplot("AT5G21150", "AGO9", ddsMat_Col)
singlegeneplot("AT3G48190", "ATM", ddsMat_Col)
singlegeneplot("AT5G40820", "ATR", ddsMat_Col)
singlegeneplot("AT1G25580", "SOG1", ddsMat_Col)

# Correlation Scatter Plot
##Relabel metadata
metadata_alt = metadata
metadata_alt[,3] = rownames(metadata)
metadata_alt[5,3] = "X68.18D.2"
metadata_alt[5,1] = "X68.18D"
metadata_alt[6,3] = "X13.6.2"
metadata_alt[6,1] = "X13.6"
metadata_alt[7,3] = "X13.6.1"
metadata_alt[7,1] = "X13.6"
metadata_alt[7,3] = "X13.6.1"
metadata_alt[7,1] = "X13.6"
metadata_alt[9,3] = "X13.9.1"
metadata_alt[9,1] = "X13.9"
metadata_alt[10,3] = "X13.9.2"
metadata_alt[10,1] = "X13.9"
metadata_alt[11,3] = "X68.21C.1"
metadata_alt[11,1] = "X68.21C"
metadata_alt[12,3] = "X68.21C.2"
metadata_alt[12,1] = "X68.21C"
metadata_alt
#Flowering time genes correlation
#read in flowering time gene data
florid_data = read.csv("/Users/emmacorcoran/Documents/R Scripts/RNA-seq/FLORID_all.csv")
#sort
florid_data = florid_data[order(florid_data$Gene.Details),]
#filter to only include flowering time genes
gene_list <- normalized_counts %>% filter(gene %in% florid_data$Gene.Details)
# Add gene labels
gene_list <- gene_list %>% mutate(genenames = florid_data$Short.Name)
# Gathering the columns to have normalized counts to a single column
gathered_list <- gene_list %>%
  gather(colnames(gene_list)[2:16], key = "sampleid", value = "normalized_counts")
#combine metadata with normalized counts
gathered_list <- inner_join(metadata_alt, gathered_list)

#Make tibbles
make_tibbles = function(results){
  toreturn_counts <- results %>% 
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  return(toreturn_counts)
}
s13.6vsCol_counts = make_tibbles(res_13.6vsCol)
s13.9vsCol_counts = make_tibbles(res_13.9vsCol)
s68.18DvsCol_counts = make_tibbles(res_68.18DvsCol)
s68.21CvsCol_counts = make_tibbles(res_68.21CvsCol)
rlt12vsCol_counts = make_tibbles(res_rlt12vsCol)
arid5vsCol_counts = make_tibbles(res_arid5vsCol)
chrvsCol_counts = make_tibbles(res_chrvsCol)

#Filter rows to only include florid flowering genes and make matrix
filter_florid = function(counts, florid_data){
  flowering_data <- counts %>% filter(gene %in% florid_data$Gene.Details)
  log_data = flowering_data$log2FoldChange
  return(log_data)
}
s13.6_flowering <- s13.6vsCol_counts %>% filter(gene %in% florid_data$Gene.Details)
s13.6_logs = filter_florid(s13.6vsCol_counts, florid_data)
s13.9_logs = filter_florid(s13.9vsCol_counts, florid_data)
s68.18D_logs = filter_florid(s68.18DvsCol_counts, florid_data)
s68.21C_logs = filter_florid(s68.21CvsCol_counts, florid_data)
rlt12_logs = filter_florid(rlt12vsCol_counts, florid_data)
arid5_logs = filter_florid(arid5vsCol_counts, florid_data)
chr_logs = filter_florid(chrvsCol_counts, florid_data)
flowering_matrix = cbind(s13.6_logs,s13.9_logs,s68.18D_logs,s68.21C_logs,rlt12_logs,arid5_logs,chr_logs)
rownames(flowering_matrix) = s13.6_flowering$gene
colnames(flowering_matrix) = c("13.6","13.9","68.18D","68.21C","rlt12","arid5","chr11chr17")
flowering_matrix = data.frame(flowering_matrix)
flowering_matrix = na.omit(flowering_matrix)

#All genes correlation
logchange_matrix = cbind(s13.6vsCol_counts$log2FoldChange,s13.9vsCol_counts$log2FoldChange,s68.18DvsCol_counts$log2FoldChange,s68.21CvsCol_counts$log2FoldChange,
                         rlt12vsCol_counts$log2FoldChange,arid5vsCol_counts$log2FoldChange,chrvsCol_counts$log2FoldChange)
rownames(logchange_matrix) = s13.6vsCol_counts$gene
colnames(logchange_matrix) = c("X13.6","X13.9","X68.18D","X68.21C","rlt12","arid5","chr11chr17")
logchange_matrix = data.frame(logchange_matrix)
logchange_matrix = na.omit(logchange_matrix)

#Make correlation scatterplot
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)
#Initiate a ggplot
b <- ggplot(logchange_matrix, aes(x = X13.6, y = X13.9, label=rownames(logchange_matrix))) + geom_point() +
  geom_smooth(method = lm, se = FALSE)
# Basic scatter plot + regression line
plotly_graph = ggplotly(b, tooltip = "all")
htmlwidgets::saveWidget(plotly_graph, "13.6vs13.9_all_correlation.html")

#Regression plot
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
#FT genes
ggplotRegression(lm(X68.18D ~ X13.6, data = flowering_matrix))
ggsave("Regression/68.18Dvs13.6_flowering_correlation.pdf")
ggplotRegression(lm(X13.9 ~ X13.6, data = flowering_matrix))
ggsave("Regression/13.9vs13.6_flowering_correlation.pdf")
ggplotRegression(lm(rlt12 ~ X13.6, data = flowering_matrix))
ggsave("Regression/rlt12vs13.6_flowering_correlation.pdf")
ggplotRegression(lm(arid5 ~ X13.6, data = flowering_matrix))
ggsave("Regression/arid5vs13.6_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ X13.6, data = flowering_matrix))
ggsave("Regression/chrvs13.6_flowering_correlation.pdf")

ggplotRegression(lm(X68.18D ~ X13.9, data = flowering_matrix))
ggsave("Regression/68.18Dvs13.9_flowering_correlation.pdf")
ggplotRegression(lm(X68.21C ~ X13.9, data = flowering_matrix))
ggsave("Regression/68.21Cvs13.9_flowering_correlation.pdf")
ggplotRegression(lm(rlt12 ~ X13.9, data = flowering_matrix))
ggsave("Regression/rlt12vs13.9_flowering_correlation.pdf")
ggplotRegression(lm(arid5 ~ X13.9, data = flowering_matrix))
ggsave("Regression/arid5vs13.9_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ X13.9, data = flowering_matrix))
ggsave("Regression/chrvs13.9_flowering_correlation.pdf")

ggplotRegression(lm(arid5 ~ rlt12, data = flowering_matrix))
ggsave("Regression/arid5vsrlt12_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ rlt12, data = flowering_matrix))
ggsave("Regression/chrvsrlt12_flowering_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ arid5, data = flowering_matrix))
ggsave("Regression/chrvsarid5_flowering_correlation.pdf")

#all genes
ggplotRegression(lm(X68.18D ~ X13.6, data = logchange_matrix))
ggsave("Regression/68.18Dvs13.6_all_correlation.pdf")
ggplotRegression(lm(X13.9 ~ X13.6, data = logchange_matrix))
ggsave("Regression/13.9vs13.6_all_correlation.pdf")
ggplotRegression(lm(rlt12 ~ X13.6, data = logchange_matrix))
ggsave("Regression/rlt12vs13.6_all_correlation.pdf")
ggplotRegression(lm(arid5 ~ X13.6, data = logchange_matrix))
ggsave("Regression/arid5vs13.6_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ X13.6, data = logchange_matrix))
ggsave("Regression/chrvs13.6_all_correlation.pdf")

ggplotRegression(lm(X68.18D ~ X13.9, data = logchange_matrix))
ggsave("Regression/68.18Dvs13.9_all_correlation.pdf")
ggplotRegression(lm(X68.21C ~ X13.9, data = logchange_matrix))
ggsave("Regression/68.21Cvs13.9_all_correlation.pdf")
ggplotRegression(lm(rlt12 ~ X13.9, data = logchange_matrix))
ggsave("Regression/rlt12vs13.9_all_correlation.pdf")
ggplotRegression(lm(arid5 ~ X13.9, data = logchange_matrix))
ggsave("Regression/arid5vs13.9_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ X13.9, data = logchange_matrix))
ggsave("Regression/chrvs13.9_all_correlation.pdf")

ggplotRegression(lm(arid5 ~ rlt12, data = logchange_matrix))
ggsave("Regression/arid5vsrlt12_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ rlt12, data = logchange_matrix))
ggsave("Regression/chrvsrlt12_all_correlation.pdf")
ggplotRegression(lm(chr11chr17 ~ arid5, data = logchange_matrix))
ggsave("Regression/chrvsarid5_all_correlation.pdf")

#####Finding Pathways for Differentially Expressed Genes
#http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

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
  ggsave(paste0("Pathways/", name, " KEGG Enrichment Pathways.pdf"), device=pdf())
  dev.off()
  #Enrich genes using Gene Ontology
  go_enrich <- enrichGO(gene = names(gene_matrix), OrgDb = 'org.At.tair.db', readable = T,
                        ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  # Plot results
  barplot(go_enrich,drop = TRUE, showCategory = 10, title = paste(name, "GO Biological Pathways"), font.size = 8)
  ggsave(paste0("Pathways/",name," GO Biological Pathways.pdf"), device=pdf())
  dev.off()
}

#Generate barplots of pathways
pathway_plots(results_sig_13.6vsCol, "13.6")
pathway_plots(results_sig_13.9vsCol, "13.9")
pathway_plots(results_sig_68.18DvsCol, "68.18D")
pathway_plots(results_sig_68.21CvsCol, "68.21C")
pathway_plots(results_sig_rlt12vsCol, "rlt12")
pathway_plots(results_sig_arid5vsCol, "arid5")
pathway_plots(results_sig_chrvsCol, "chr11chr17")
