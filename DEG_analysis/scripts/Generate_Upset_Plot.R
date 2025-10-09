#### Libraries #####
library(tidyverse)
library("UpSetR")
library(ComplexUpset)
library(devtools)
library(ggforce)

#### Functions ####

load_matrix <- function(filename){
  return(read_csv(file = filename, col_names = TRUE,))
}

save_plots <- function(plot_obj, width = 8, height = 6, graph_name = "Rplot"){
  eps_name = paste0(graph_name, ".eps")
  pdf_name = paste0(graph_name, ".pdf")
  
  setEPS(width = width, height = height)
  postscript(eps_name)
  print(plot_obj)
  dev.off()
  
  pdf(file = pdf_name, width = width, height = 6)
  print(plot_obj)
  dev.off()
}


get_transcripts <- function(dir_path, cond){
  DGE_table = read_tsv(paste0(dir_path, cond, "/DGE_result_",cond,"_FC_and_FDR_filter.tsv"))
  return(DGE_table[,1])
}

get_DGE_from_list <- function(dir_path, cond, list_of_transcript){
  DGE_table = read_tsv(paste0(dir_path, cond, "/DGE_result_FC_and_FDR_filter_fully_annotated.tsv"))
  res = DGE_table[which(DGE_table$X.. %in% list_of_transcript),]
  dir_path = "D:/Louis/Transcriptome/Regeneration_analysis/analysis/DEG_sequential/"
  write_tsv(x = res,file = paste0(dir_path, cond, "/UNIQUE_SEQUENTIEL",cond,".tsv"))
}
#### Main ####


# Conditions by conditions ----

# Loading data ----

matrix_NA_1TMM = load_matrix("Table/NA/NA_1TMM_threshold_matrix.csv")

matrix_0dpa_1TMM = load_matrix("Table/0dpa/0dpa_1TMM_threshold_matrix.csv")

matrix_1dpa_1TMM = load_matrix("Table/1dpa/1dpa_1TMM_threshold_matrix.csv")

matrix_3dpa_1TMM = load_matrix("Table/3dpa/3dpa_1TMM_threshold_matrix.csv")

matrix_2dpa_1TMM = load_matrix("Table/2dpa/2dpa_1TMM_threshold_matrix.csv")

matrix_5dpa_1TMM = load_matrix("Table/5dpa/5dpa_1TMM_threshold_matrix.csv")


# Upset plot ----

y <- list(
  "0 dpa" = as.data.frame(matrix_0dpa_1TMM)[,1],
  "1 dpa" = as.data.frame(matrix_1dpa_1TMM)[,1],
  "2 dpa" = as.data.frame(matrix_3dpa_1TMM)[,1],
  "3 dpa" = as.data.frame(matrix_2dpa_1TMM)[,1],
  "5 dpa" = as.data.frame(matrix_5dpa_1TMM)[,1],
  "NA" = as.data.frame(matrix_NA_1TMM)[,1]
)


Upset_plot = UpSetR::upset(fromList(y),number.angles = 0, group.by = "degree", order.by = "degree", decreasing = F,
      line.size = 0.5,point.size = 2, nintersects = NA,
      mainbar.y.label = "Gene intersections", sets.x.label = "Genes Per Samples",
      sets = c("NA","5 dpa","3 dpa","2 dpa","1 dpa","0 dpa"),
      text.scale = c(1.3,1.3,1,1,2,1.1), keep.order = T, nsets = 6,
      matrix.color = "black", main.bar.color = "black", show.numbers = "yes",
      sets.bar.color =  "grey", mainbar.y.max = 20000, #scale.intersections = "log2",
      queries = list(list(query= intersects, params = list("NA"), color = "orange", active = T),
                     list(query=intersects, params = list("0 dpa"), color = "orange", active = T),
                     list(query=intersects, params = list("1 dpa"), color = "orange", active = T),
                     list(query=intersects, params = list("2 dpa"), color = "orange", active = T),
                     list(query=intersects, params = list("3 dpa"), color = "orange", active = T),
                     list(query=intersects, params = list("5 dpa"), color = "orange", active = T)))


save_plots(Upset_plot, width = 10, height = 7,
           graph_name = "Figures/Upset_plot_1TMM_limit")


# Based on sequential DGE. ----

conditions = c("0dpavs1dpa", "1dpavs2dpa", "2dpavs3dpa", "3dpavs5dpa", "5dpavsNA")

dir_path = "DEG/"

y <- list()

for (cond in conditions) {
  y[cond] = get_transcripts(dir_path, cond)
}

Upset_plot = UpSetR::upset(fromList(y),number.angles = 0, group.by = "degree", order.by = "degree", decreasing = F,
                           line.size = 0.5,point.size = 2, nintersects = NA,
                           mainbar.y.label = "Gene intersections", sets.x.label = "Genes per conditions",
                           sets = rev(conditions),
                           text.scale = c(1.3,1.3,1,1,2,1.1), keep.order = T, nsets = 5,
                           matrix.color = "black", main.bar.color = "black", show.numbers = "yes",
                           sets.bar.color =  "grey", mainbar.y.max = 4500, #scale.intersections = "log2",
                           queries = list(list(query= intersects, params = list("0dpavs1dpa"), color = "orange", active = T),
                                          list(query=intersects, params = list("1dpavs2dpa"), color = "orange", active = T),
                                          list(query=intersects, params = list("2dpavs3dpa"), color = "orange", active = T),
                                          list(query=intersects, params = list("3dpavs5dpa"), color = "orange", active = T),
                                          list(query=intersects, params = list("5dpavsNA"), color = "orange", active = T)))

Upset_plot

save_plots(Upset_plot, width = 10, height = 7,
           graph_name = "figures/Upset_plot_DGE_freq")


# Save lists corresponding to unique sequential conditions ####

dpa0vsdpa1 = setdiff(x = y[[1]], y = c(y[[2]],y[[3]],y[[4]],y[[5]]))
dpa1vsdpa2 = setdiff(x = y[[2]], y = c(y[[1]],y[[3]],y[[4]],y[[5]]))
dpa2vsdpa3 = setdiff(x = y[[3]], y = c(y[[1]],y[[2]],y[[4]],y[[5]]))
dpa3vsdpa5 = setdiff(x = y[[4]], y = c(y[[1]],y[[2]],y[[3]],y[[5]]))
dpa5vsNA =  setdiff(x = y[[5]], y = c(y[[1]],y[[2]],y[[3]],y[[4]]))

get_DGE_from_list(dir_path = dir_path, cond = conditions[1], list_of_transcript = dpa0vsdpa1)
get_DGE_from_list(dir_path = dir_path, cond = conditions[2], list_of_transcript = dpa1vsdpa2)
get_DGE_from_list(dir_path = dir_path, cond = conditions[3], list_of_transcript = dpa2vsdpa3)
get_DGE_from_list(dir_path = dir_path, cond = conditions[4], list_of_transcript = dpa3vsdpa5)
get_DGE_from_list(dir_path = dir_path, cond = conditions[5], list_of_transcript = dpa5vsNA)




####


list_of_genes = matrix_NA_1TMM$...1[which( ! matrix_NA_1TMM$...1 %in% matrix_0dpa_1TMM$...1 &
              ! matrix_NA_1TMM$...1 %in% matrix_1dpa_1TMM$...1 &
              ! matrix_NA_1TMM$...1 %in% matrix_2dpa_1TMM$...1 &
              ! matrix_NA_1TMM$...1 %in% matrix_3dpa_1TMM$...1 &
              ! matrix_NA_1TMM$...1 %in% matrix_5dpa_1TMM$...1)]
write.csv(x = list_of_genes, file = "figures/genes_uniques/NA.csv",
            quote = F, row.names = F)

conds = c("0dpa", "1dpa", "2dpa", "3dpa", "5dpa", "NA")

mouse_annotation = read_tsv(
  file = "Mus_musculus_table.tsv",
  col_names = TRUE)

for (i in conds) {
  actual_cond = read_csv(file = paste0("figures/genes_uniques/",i,".csv"),
                         col_names = F)
  blast_name = c()
  for (y in actual_cond$X1) {
    blast_name = c(blast_name, mouse_annotation$Name[which(mouse_annotation$TranscriptID == paste0(y, "t1"))])
  }
  actual_cond$mouse = blast_name
  write_csv(x = actual_cond, file = paste0("figures/genes_uniques/", i, "_annotation.csv"))
}
