# Code to execute DGE analysis on the regeneration samples (single-end).
# Several steps are done :
# Filtering, quality control and then statistical anaylsis for DGE analysis and plotting (volcano plot ..)


# Libraries ###
library(edgeR)
library(rgl)
library(ggfortify)
library("FactoMineR")
library("factoextra")
library("smacof")
library("dplyr")
library(ggplot2)
library(gplots)
library(statmod)
library(magrittr)
library(EnhancedVolcano)
# Seurat library for dark mode.
library(Seurat)
library(tidyverse)
library("dendextend")

# Functions ####

# Plot a 3D pca.
plotPCA <- function(x, groups) {
  plot3d(x, col=c("red","blue","orange","green","yellow","grey"), type="s", size=1, axes=F)
  axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
  grid3d("x")
  grid3d("y")
  grid3d("z")
}

# Plot a sample correlation matrix using the log CPM count and a list
# giving column indexes.
plot_individual_cor_heatmap <- function(matrix, list_pos){
  
  nbr_of_rep = length(list_pos)
  head_names = colnames(matrix[,list_pos])
  submatrix = matrix[,list_pos]
  cor_matrix = matrix(nrow = nbr_of_rep, ncol = nbr_of_rep)
  rownames(cor_matrix) = head_names
  colnames(cor_matrix) = head_names
  
  for (i in 1:nbr_of_rep) {
    for (y in 1:nbr_of_rep) {
      cor_matrix[i,y] = cor(x = submatrix[,i], y = submatrix[,y])
    }
  }
  colors = c(seq(0.5,0.7,length=100),seq(0.71,0.79,length=100), seq(0.80,1,length=100))
  mypalette = colorRampPalette(c("purple","black","yellow"))(n = 299)
  
  par(mar = c(5.1, 4.1, 4.1, 10))
  heatmap.2(cor_matrix, col=mypalette,breaks = colors, key=T, keysize=1.5,
            density.info="none", trace="none",cexCol=0.9, cexRow = 0.9,
            main = "Sample correlation matrix", margins = c(6.5,7),
            dendrogram = 'both', symkey = F, symm = F, symbreaks = F,
            scale = "none", lwid = c(1,5), lhei = c(1.5,5)
  )
  
}


# Shortened the output from trinotate annotation by taking only the first
# result from a line.
shorter <- function(str2short){
  if(grepl("^", str2short, fixed = T)){
    return(unlist(strsplit(str2short,split = "^", fixed = T))[[1]])
  }else{
    return(unlist(strsplit(str2short,split = "`", fixed = T))[[1]])
  }
}


# Concatenate a list into a single string.
agglom <- function(str2aggl){
  aggl = ""
  for (y in 1:length(unlist(str2aggl))) {
    aggl = paste(aggl, shorter(unlist(str2aggl)[y]))
  }
  return(aggl)
}

# Add annotation to a Edger DGE or DTE output matrix using the trinotate report.
add_annotation <- function(table, annotation){
  
  blastx = c()
  blastp = c()
  GOterms = c()
  kegg = c()
  
  for (i in 1:length(table[,1])) {
    line = annotation[which(annotation$X.gene_id == rownames(table)[i]),]
    blastx = c(blastx, agglom(line[3]))
    blastp = c(blastp, agglom(line[7]))
    GOterms = c(GOterms, agglom(line[13]))
    kegg = c(kegg, agglom(line[12])) 
  }
  table$blastx = blastx
  table$blastp = blastp
  table$GOterms = GOterms
  table$kegg = kegg
  return(table)
}

##


# Main ####

# Data Loading ----
# results from trinotate
trinotate = read.table("trinotate_annotation_report.xls",
                       header=T, stringsAsFactors = F,sep = "\t", quote="", com='')

# raw count from kallisto
data = read.table("data/kallisto.gene.counts.matrix",
                  header = T, row.names = 1, com='', sep = "\t")

# Contrast preparation ----

colnames(data)

# remove outlier, ignore this 3 lines to see analysis with outliers (need to adapt code afterhand)
data = data[,-18] # 5dpa rep3
data = data[,-13] # 3dpa rep1
data = data[,-9] # 1dpa rep3

data = data[,c(4:15,1,2,3)]

rep = factor(c("rep1","rep2","rep3",
               "rep1","rep2",
               "rep1","rep2","rep3",
               "rep2","rep3",
               "rep1","rep2",
               "rep1","rep2","rep3"))

cond = factor(c("0dpa", "0dpa", "0dpa",
                "1dpa", "1dpa",
                "2dpa","2dpa", "2dpa",
                "3dpa","3dpa",
                "5dpa","5dpa",
                "NA", "NA", "NA"))

colors_conds_no_filter = c("#66C2A5","#66C2A5","#66C2A5",
  "#FC8D62","#FC8D62","#FC8D62",
                                    "#8DA0CB","#8DA0CB","#8DA0CB",
                                    "#E78AC3","#E78AC3","#E78AC3",
                                    "#A6D854","#A6D854","#A6D854",
           "#FFD92F","#FFD92F","#FFD92F")
                                    
colors_conds = c("#66C2A5","#66C2A5","#66C2A5",
                          "#FC8D62","#FC8D62",
                          "#8DA0CB","#8DA0CB","#8DA0CB",
                          "#E78AC3","#E78AC3",
                          "#A6D854","#A6D854",
                          "#FFD92F","#FFD92F","#FFD92F")

cond = relevel(cond, ref="0dpa")


# Create EdgeR object and filter out lowly or non expressed genes ----

matrix_NA_1TMM = read_csv(file = "Table/NA/NA_1TMM_threshold_matrix.csv")
matrix_0dpa_1TMM = read_csv(file = "Table/0dpa/0dpa_1TMM_threshold_matrix.csv")
matrix_1dpa_1TMM = read_csv(file = "Table/1dpa/1dpa_1TMM_threshold_matrix.csv")
matrix_2dpa_1TMM = read_csv(file = "Table/2dpa/2dpa_1TMM_threshold_matrix.csv")
matrix_3dpa_1TMM = read_csv(file = "Table/3dpa/3dpa_1TMM_threshold_matrix.csv")
matrix_5dpa_1TMM = read_csv(file = "Table/5dpa/5dpa_1TMM_threshold_matrix.csv")

keep = c(matrix_NA_1TMM$...1,
         matrix_0dpa_1TMM$...1,
         matrix_1dpa_1TMM$...1,
         matrix_2dpa_1TMM$...1,
         matrix_3dpa_1TMM$...1,
         matrix_5dpa_1TMM$...1) %>% unique()

all_gene = rownames(data)
keep2 = rep(FALSE,length(all_gene))
names(keep2) = all_gene
for (i in 1:length(all_gene)) {
  if(all_gene[i] %in% keep){
    keep2[i] = TRUE
  }
}

data = DGEList(counts = data, group = cond)
keep = filterByExpr(data)
table(keep2)
data <- data[keep2, , keep.lib.sizes = F]


data = calcNormFactors(data)

### Analysis ----

## Data quality ----

# Batch effect?
# Plot MDS

data_mds = edgeR::cpm(data, log=TRUE, prior.count = 1)
plotMDS(data_mds, col=c(1,1,1,2,2,2,3,3,4,4,4,5,5,6,6), main = "MDS plot before batch correction")
plotMDS(data_mds, col=colors_conds, main = "MDS plot before batch correction")

# remove batch effect
data_bc = removeBatchEffect(cpm(data, log=TRUE, prior.count = 1), batch = rep)
plotMDS(data_bc, col=rep(1:2, each = 3), main = "MDS plot after batch correction")

# Batch effect correction not needed.


# Compute top 1000 variable gene

data_sd = apply(data_mds, MARGIN = 1, FUN = sd)

data_sd = data.frame(data_mds, data_sd)

data_sd = data_sd[order(-data_sd$data_sd),]

top1000 = head(data_sd, n = 1000)

# PCA generation ----
# 2D PCA
res.pca <- PCA(t(data_mds), graph = F)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,100))
var <- get_pca_var(res.pca)
fviz_pca_ind(res.pca, geom = c("point","text"),
             col.var = c("0dpa", "0dpa", "0dpa",
                         "1dpa", "1dpa",
                         "2dpa","2dpa", "2dpa",
                         "3dpa", "3dpa",
                         "5dpa", "5dpa",
                         "NA", "NA", "NA"),
             title = "", col.ind = cond, addEllipses = F,
             palette = colors_conds_no_filter[seq(from = 1,
                                                  to = length(colors_conds_no_filter),
                                                  by = 3)]) + 
  theme(panel.background = element_rect(fill= "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_shape_manual(values=c(15,15,15,15,15,15))

# 3D PCA -> labeling problem..
pca = prcomp(t(data_mds), scale. = T)
plotPCA(pca$x[,1:3], cond)

# Dendrogram ----


dd <- dist(t(data_mds), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

fit.dendrogram <- as.dendrogram(hc)

labels_colors(fit.dendrogram) <- colors_conds[order.dendrogram(fit.dendrogram)]

plot(fit.dendrogram, main = "Dendrogram clustering", xlab = "")


# Heatmap generation ----

# Basic heatmap -> using top 1000, else too much genes
heatmap(x = as.matrix(top1000), main = "Heatmap",
        ylab = "Genes", xlab = "Condition")

# Heatmap with heatmap2
# matrix arrangement for heatmap.
mat_df =as.matrix(top1000[,1:15])
rownames(mat_df) = rownames(top1000)

colors = c(seq(-10,-2,length=100),seq(-1.9,1.9,length=100), seq(2,10,length=100))
mypalette = colorRampPalette(c("white","orange","red"))(n = 299)


heatmap.2(mat_df, col=mypalette,breaks = colors, key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9,
          main = "Heatmap of the top 1000 variable transcripts", labRow = F, margins = c(6.5,1.5),
          dendrogram = 'column', symkey = F, symm = F, symbreaks = F,
          scale = "none", xlab = "Samples", ylab = "Transcripts")

## Correlation matrix ----

cor_matrix = matrix(nrow = length(rep), ncol = length(rep))
rownames(cor_matrix) = paste(cond,rep,sep="_")
colnames(cor_matrix) = paste(cond,rep,sep="_")
cor_matrix
for (i in 1:length(cor_matrix[,1])) {
  for (y in 1:length(cor_matrix[,1])) {
    cor_matrix[i,y] = cor(x = data_mds[,i], y = data_mds[,y])
  }
}

colors = c(seq(0.5,0.7,length=10),seq(0.71,0.79,length=10), seq(0.80,1,length=10))
mypalette = colorRampPalette(c("purple","black","yellow"))(n = 29)

par(mar = c(5.1, 4.1, 4.1, 3.1))
heatmap.2(cor_matrix, col=mypalette,breaks = colors, key=T, keysize=1.5,
          density.info="none", trace="none",cexCol=0.9, Rowv = T, Colv = T,
          main = "Sample correlation matrix", margins = c(6.5,7),
          dendrogram = 'both', symkey = F, symm = F, symbreaks = F,
          scale = "none", lwid = c(1,5), lhei = c(1.5,5),
          RowSideColors = c(
            colors_conds
          ),
          ColSideColors = c(
            colors_conds
          )
        )



###### Differential expression ----

# contrast
design = model.matrix(~0+cond)
rownames(design) <- colnames(data)
data <- estimateDisp(data, design, robust = T)

######### Automatic pairwise tests ----

cond_fixed = unique(paste("cond", cond, sep =""))
nbr_of_loop = 0
nbr_of_cond = length(cond_fixed)
for (i in 1:(nbr_of_cond-1)){
  for(y in 1:(nbr_of_cond-i)){
    ## Contrast making
    nbr_of_loop = nbr_of_loop+1
    cond1 = cond_fixed[i]
    cond2 = cond_fixed[i+y]
    actual_contrast = makeContrasts(paste(cond2, "-", cond1,sep=""), levels = cond_fixed)
    ##
    
    if (!dir.exists(paste0("DEG/",substr(cond1, 5, 100),"vs",substr(cond2, 5, 100)))){
    dir.create(paste0("DEG/",substr(cond1, 5, 100),"vs",substr(cond2, 5, 100)))
    }
    
    print(actual_contrast)
    
    ## Statistics tests :
    lrt <- glmLRT(fit, contrast = actual_contrast)
    DET_matrix = topTags(lrt, n = length(lrt$table[,1]))
    ##
    
    ## Save volcano plot
    volcano = EnhancedVolcano(DET_matrix$table,
                              lab = NA,
                              x = 'logFC',
                              y = 'FDR',
                              #ylim = c(0,20),
                              title = paste("Volcano Plot of ",cond1, " vs ", cond2, " conditions",sep=""),
                              subtitle = bquote(italic("Log2 FC threshold = 1 and p_value (FDR) threshold = 0.05")),
                              FCcutoff = 1,
                              pCutoff = 0.05,
                              caption = paste0("total = ", nrow(DET_matrix$table), " genes"))
    ggsave(filename = paste("DEG/",substr(cond1, 5, 100),"vs",substr(cond2, 5, 100),"/volcano_plot_test.pdf", sep = ""),
           plot = volcano, device = "pdf", height = 6.5,width = 7)
    
    ggsave(filename = paste("DEG/",substr(cond1, 5, 100),"vs",substr(cond2, 5, 100),"/volcano_plot_test.png", sep = ""),
           plot = volcano, device = "png", height = 6.5,width = 7)
    ##
    
    ## Apply thresholds
    thresholded_matrix = DET_matrix$table[which(abs(DET_matrix$table$logFC) > 1 & DET_matrix$table$FDR < 0.05),]
    
    ## Save MA plot
    pdf(file = paste("DEG/",substr(cond1, 5, 100),"vs",substr(cond2, 5, 100),"/MA_plot_test.pdf", sep = ""))
    plotSmear(lrt, de.tags = rownames(thresholded_matrix))
    dev.off()
    ##
    
    ## Save table with threshold
    write.table(x = thresholded_matrix,
                file = paste("DEG/",substr(cond1, 5, 100),"vs",
                             substr(cond2, 5, 100),"/DGE_result_",
                             substr(cond1, 5, 100),"vs",
                             substr(cond2, 5, 100),"_FC_and_FDR_filter.tsv", sep = ""),
                sep = "\t", row.names = T, col.names = NA)
    ## Save table without threshold
    write.table(x = DET_matrix$table, 
                file = paste("DEG/",substr(cond1, 5, 100),
                             "vs",substr(cond2, 5, 100),
                             "/DGE_result_",substr(cond1, 5, 100),
                             "vs",substr(cond2, 5, 100),".tsv", sep = ""),
                sep = "\t", row.names = T, col.names = NA)
    
    annotated_matrix = add_annotation(thresholded_matrix, trinotate)
    
    write.table(x = annotated_matrix, file = paste("DEG/",substr(cond1, 5, 100),"vs",
                                                   substr(cond2, 5, 100),"/DGE_result_",
                                                   substr(cond1, 5, 100),"vs",
                                                   substr(cond2, 5, 100),"_FC_and_FDR_filter_trinotate.tsv", sep=""),
                sep = "\t", row.names = T, col.names = NA)
    
  }
}


###########################