# Libraries
library("VennDiagram")
library("UpSetR")
library(ComplexUpset)
library(ggvenn)
library(ggVennDiagram)
library(ggplot2)

#######
# /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
# /!\ Input filename (TPM_table) has been replaced by the TMM matrix /!\
# /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
#######

### Functions

save_matrix <- function(matrix, filename){
  write.table(x = matrix, file = filename, sep = ",", row.names = T, col.names = NA)
}

####

### Main

# Loading data
TMM_table = read.table("data/kallisto.gene.TMM.EXPR.matrix",
                       header=T, row.names=1, com='', check.names=F)

trinotate = read.table("data/trinotate_annotation_report.xls",
                       header=T, stringsAsFactors = F,sep = "\t", quote="", com='')
#

# Add median information
TMM_table$median_NA = apply(TMM_table[,c(1,2,3)],1,median)
TMM_table$median_0dpa = apply(TMM_table[,c(4,5,6)],1,median)
TMM_table$median_1dpa = apply(TMM_table[,c(7,8)],1,median)
TMM_table$median_2dpa = apply(TMM_table[,c(10,11,12)],1,median)
TMM_table$median_3dpa = apply(TMM_table[,c(14,15)],1,median)
TMM_table$median_5dpa = apply(TMM_table[,c(16,17)],1,median)

#

## NA matrix generation
# All genes
matrix_NA = TMM_table[,c(1,2,3,19)]
save_matrix(matrix_NA, "Table/NA/NA_TMM_matrix.csv")

# 1TMM threshold
matrix_NA_1TMM = matrix_NA[which(matrix_NA$median_NA>1),]
save_matrix(matrix_NA_1TMM, "Table/NA/NA_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_NA_5TMM = matrix_NA[which(matrix_NA$median_NA>5),]
save_matrix(matrix_NA_5TMM, "Table/NA/NA_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_NA_10TMM = matrix_NA[which(matrix_NA$median_NA>10),]
save_matrix(matrix_NA_10TMM, "Table/NA/NA_10TMM_threshold_matrix.csv")

#

## 0dpa matrix generation
# All genes
matrix_0dpa = TMM_table[,c(4,5,6,20)]
save_matrix(matrix_0dpa, "Table/0dpa/0dpa_TMM_matrix.csv")

# 1TMM threshold
matrix_0dpa_1TMM = matrix_0dpa[which(matrix_0dpa$median_0dpa>1),]
save_matrix(matrix_0dpa_1TMM, "Table/0dpa/0dpa_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_0dpa_5TMM = matrix_0dpa[which(matrix_0dpa$median_0dpa>5),]
save_matrix(matrix_0dpa_5TMM, "Table/0dpa/0dpa_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_0dpa_10TMM = matrix_0dpa[which(matrix_0dpa$median_0dpa>10),]
save_matrix(matrix_0dpa_10TMM, "Table/0dpa/0dpa_10TMM_threshold_matrix.csv")

#

## 1dpa matrix generation
# All genes
matrix_1dpa = TMM_table[,c(7,8,21)]
save_matrix(matrix_1dpa, "Table/1dpa/1dpa_TMM_matrix.csv")

# 1TMM threshold
matrix_1dpa_1TMM = matrix_1dpa[which(matrix_1dpa$median_1dpa>1),]
save_matrix(matrix_1dpa_1TMM, "Table/1dpa/1dpa_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_1dpa_5TMM = matrix_1dpa[which(matrix_1dpa$median_1dpa>5),]
save_matrix(matrix_1dpa_5TMM, "Table/1dpa/1dpa_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_1dpa_10TMM = matrix_1dpa[which(matrix_1dpa$median_1dpa>10),]
save_matrix(matrix_1dpa_10TMM, "Table/1dpa/1dpa_10TMM_threshold_matrix.csv")

#


## 2dpa matrix generation
# All genes
matrix_2dpa = TMM_table[,c(10,11,12,22)]
save_matrix(matrix_2dpa, "Table/2dpa/2dpa_TMM_matrix.csv")

# 1TMM threshold
matrix_2dpa_1TMM = matrix_2dpa[which(matrix_2dpa$median_2dpa>1),]
save_matrix(matrix_2dpa_1TMM, "Table/2dpa/2dpa_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_2dpa_5TMM = matrix_2dpa[which(matrix_2dpa$median_2dpa>5),]
save_matrix(matrix_2dpa_5TMM, "Table/2dpa/2dpa_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_2dpa_10TMM = matrix_2dpa[which(matrix_2dpa$median_2dpa>10),]
save_matrix(matrix_2dpa_10TMM, "Table/2dpa/2dpa_10TMM_threshold_matrix.csv")

#

## 3dpa matrix generation
# All genes
matrix_3dpa = TMM_table[,c(14,15,23)]
save_matrix(matrix_3dpa, "Table/3dpa/3dpa_TMM_matrix.csv")

# 1TMM threshold
matrix_3dpa_1TMM = matrix_3dpa[which(matrix_3dpa$median_3dpa>1),]
save_matrix(matrix_3dpa_1TMM, "Table/3dpa/3dpa_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_3dpa_5TMM = matrix_3dpa[which(matrix_3dpa$median_3dpa>5),]
save_matrix(matrix_3dpa_5TMM, "Table/3dpa/3dpa_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_3dpa_10TMM = matrix_3dpa[which(matrix_3dpa$median_3dpa>10),]
save_matrix(matrix_3dpa_10TMM, "Table/3dpa/3dpa_10TMM_threshold_matrix.csv")



## 5dpa matrix generation
# All genes
matrix_5dpa = TMM_table[,c(16,17,24)]
save_matrix(matrix_5dpa, "Table/5dpa/5dpa_TMM_matrix.csv")

# 1TMM threshold
matrix_5dpa_1TMM = matrix_5dpa[which(matrix_5dpa$median_5dpa>1),]
save_matrix(matrix_5dpa_1TMM, "Table/5dpa/5dpa_1TMM_threshold_matrix.csv")

# 5TMM threshold
matrix_5dpa_5TMM = matrix_5dpa[which(matrix_5dpa$median_5dpa>5),]
save_matrix(matrix_5dpa_5TMM, "Table/5dpa/5dpa_5TMM_threshold_matrix.csv")

# 10TMM threshold
matrix_5dpa_10TMM = matrix_5dpa[which(matrix_5dpa$median_5dpa>10),]
save_matrix(matrix_5dpa_10TMM, "Table/5dpa/5dpa_10TMM_threshold_matrix.csv")

#
#### Upset plot

y <- list(
  dpa0 = rownames(matrix_0dpa_1TMM),
  dpa1 = rownames(matrix_1dpa_1TMM),
  dpa2 = rownames(matrix_2dpa_1TMM),
  dpa3 = rownames(matrix_3dpa_1TMM),
  dpa5 = rownames(matrix_5dpa_1TMM),
  NonA = rownames(matrix_NA_1TMM)
)

UpSetR::upset(fromList(y), number.angles = 0,
      line.size = 1.5,point.size = 3.5, order.by = c("degree"), decreasing = F,
      mainbar.y.label = "Gene intersections", sets.x.label = "Genes Per Samples",
      sets = c("NonA","dpa5","dpa3","dpa2","dpa1","dpa0"), nintersects = 70,
      text.scale = c(1.3,1.3,1,1,2,1.1), keep.order = T, nsets = 6,
      matrix.color = "black", main.bar.color = "black", show.numbers = "no",
      queries = list(list(query=intersects, params = list("dpa0"), color = "orange", active = T),
                     list(query=intersects, params = list("dpa1"), color = "orange", active = T),
                     list(query=intersects, params = list("dpa2"), color = "orange", active = T),
                     list(query=intersects, params = list("dpa3"), color = "orange", active = T),
                     list(query=intersects, params = list("dpa5"), color = "orange", active = T),
                     list(query= intersects, params = list("NonA"), color = "orange", active = T)))




####