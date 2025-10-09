# Hierarchical clustering of gene expression
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html

# Libraries ####
library(tidyverse)
library(RColorBrewer)

# Parameters ####

nbr_of_cluster = 15

# Function ####

remove_na = function(m){
  m.nona = m[!rowSums(!is.finite(m)),]
  return(m.nona)
}

save_plots <- function(plot_obj, width = 8, height = 6, graph_name = "Rplot"){
  pdf_name = paste0(graph_name, ".pdf")

  pdf(file = pdf_name, width = width, height = height)
  print(plot_obj)
  dev.off()
}


# Main ####

TMM_matrix = read.table("data/kallisto.gene.TMM.EXPR.matrix",
                        header = T, row.names = 1, com='', sep = "\t")

# remove outliers
TMM_matrix = TMM_matrix[,-18] # 5dpa rep3
TMM_matrix = TMM_matrix[,-13] # 3dpa rep1
TMM_matrix = TMM_matrix[,-9] # 1dpa rep3

median_TMM_matrix = data.frame(dpa0 = apply(X = TMM_matrix[,4:6], MARGIN = 1, FUN = median),
                               dpa1 = apply(X = TMM_matrix[,7:8], MARGIN = 1, FUN = median),
                               dpa2 = apply(X = TMM_matrix[,9:11], MARGIN = 1, FUN = median),
                               dpa3 = apply(X = TMM_matrix[,12:13], MARGIN = 1, FUN = median),
                               dpa5 = apply(X = TMM_matrix[,14:15], MARGIN = 1, FUN = median),
                               NonA = apply(X = TMM_matrix[,1:3], MARGIN = 1, FUN = median))


cond = c("0dpavs1dpa","0dpavs2dpa","0dpavs3dpa","0dpavs5dpa","0dpavsNA",
         "1dpavs2dpa","1dpavs3dpa","1dpavs5dpa","1dpavsNA",
         "2dpavs3dpa","2dpavs5dpa","2dpavsNA",
         "3dpavs5dpa","3dpavsNA",
         "5dpavsNA")

genes_id = c()

file_path = "DEG/"

for (i in cond) {
  sign_genes = read.table(paste(file_path,i,"/DGE_result_",i,"_FC_and_FDR_filter.tsv",sep = ""),
                          header = T, row.names = 1, com='', sep="\t")
  genes_id = c(genes_id,rownames(sign_genes))
}

unique_genes_id = unique(genes_id)

matrix_sign_genes = median_TMM_matrix[which(rownames(median_TMM_matrix) %in% unique_genes_id),]

#matrix_sign_genes = matrix_sign_genes[-c(1:100),]

tail(matrix_sign_genes)

# Scale matrix (gene names have to be in columns)

# ncol = 6
scaled_matrix = matrix_sign_genes %>% t() %>% scale() %>% t() %>%
  as.matrix(nrow =nrow(matrix_sign_genes), ncol = ncol(matrix_sign_genes)) %>%
  as.data.frame()

# Calculate distance

scaled_matrix = as.data.frame(scaled_matrix)
gene_dist <- dist(scaled_matrix)

# Hierarchical clustering

gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = F)
# h = 4.2
abline(h = 3.7, col = "brown", lwd = 2)


for (i in 1:nbr_of_cluster) {
  
  
  dir.path = paste("clustering/hclust_k",as.character(i),sep="")
  if(!dir.exists(dir.path)){
    dir.create(dir.path) 
  }
  
  gene_cluster = cutree(gene_hclust, k= i) %>% enframe()
  colnames(gene_cluster) = c("gene","cluster")
  ######## Visualization
  
  scaled_matrix$cluster = as.data.frame(gene_cluster[which(gene_cluster$gene == rownames(scaled_matrix)),2])
  
  ##rename clusters
  
  new_clusters_name = c()
  for (y in 1:i) {
    new_clusters_name = c(new_clusters_name, paste("Cluster n ",
                                                   as.character(y),
                                                   " : ",
                                                   as.character(length(which(scaled_matrix$cluster == y))),
                                                   " genes",
                                                   sep = ""))
  }
  
  rename_clusters = function(cluster_number, list_new_cluster_name){
    return(list_new_cluster_name[cluster_number])
  }
  
  scaled_matrix$cluster = apply(X = scaled_matrix$cluster, MARGIN = 1, FUN = rename_clusters, new_clusters_name)
  
  scaled_matrix$cluster = factor(scaled_matrix$cluster, levels = new_clusters_name)
  
  #
  
  scaled_matrix$geneID = rownames(scaled_matrix)
  
  tribble_scaled_matrix = scaled_matrix %>% pivot_longer(cols = c(1,2,3,4,5,6),names_to = "Stade")
  
  tribble_scaled_matrix$Stade = as.factor(tribble_scaled_matrix$Stade)
  
  tribble_scaled_matrix$Stade = ordered(tribble_scaled_matrix$Stade, levels = c("dpa0", "dpa1", "dpa2", "dpa3", "dpa5", "NonA"))
  
  df = as.data.frame(tribble_scaled_matrix)
  
  cluster_plot = tribble_scaled_matrix %>% ggplot(aes(Stade, value)) +
    geom_line(aes(group = geneID), alpha = 0.3, colour = "grey") +
    geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5,
              aes(group = 1))+
    facet_wrap(facets = vars(cluster), nrow = i%/%5+i%%5, ncol = 5) + theme_classic()
  ggsave(paste(dir.path, "/cluster_profile_ordered.pdf", sep = ""),
         plot = cluster_plot, width = 15, height = 12, units = "in")
  ggsave(paste(dir.path, "/cluster_profile_ordered.png", sep = ""),device = "png",
         plot = cluster_plot, width = 15, height = 12, units = "in")
  
  
  ## Save geneID for each clusters.
  
  clusters = 1:i
  scaled_matrix$cluster = as.data.frame(gene_cluster[which(gene_cluster$gene == rownames(scaled_matrix)),2])
  cluster_numbers = c()
  genes_name = c()
  for (y in clusters) {
    geneIDs = rownames(scaled_matrix[which(scaled_matrix$cluster == y),])
    cluster_numbers = c(cluster_numbers, rep(y, times = length(geneIDs)))
    genes_name = c(genes_name, geneIDs)
    name = ifelse(y < 10, paste("0",as.character(y),sep=""),y)
    dir.path= paste("analysis/clustering2/hclust_k",as.character(i),"/cluster",name,sep="")
    if(!dir.exists(dir.path)){
      dir.create(dir.path) 
    }
    write.table(x = geneIDs,
                file = paste(dir.path,"/cluster_gene_list.tsv", sep = ""),
                sep = "\t", row.names = F, col.names = F, quote = F)
  }
  df = data.frame(gene = genes_name, cluster = cluster_numbers)
  write.table(x = df,
              file = paste0("analysis/clustering2/hclust_k",as.character(i),"/gene_to_cluster_correspondance.tsv"),
              quote = F, sep = "\t", col.names = T, row.names = F)
}

# add quantile to clusters profiles

ymin_1 = c()
ymin_2 = c()
ymin_3 = c()
ymin_4 = c()

ymax_1 = c()
ymax_2 = c()
ymax_3 = c()
ymax_4 = c()

for (i in 1:nrow(tribble_scaled_matrix)) {
  local_mean = median(tribble_scaled_matrix$value[which(tribble_scaled_matrix$cluster == tribble_scaled_matrix$cluster[i] &
                                                          tribble_scaled_matrix$Stade == tribble_scaled_matrix$Stade[i])])
  
  quantile_1 = quantile(tribble_scaled_matrix$value[which(tribble_scaled_matrix$cluster == tribble_scaled_matrix$cluster[i] &
                                                            tribble_scaled_matrix$Stade == tribble_scaled_matrix$Stade[i] & 
                                                            tribble_scaled_matrix$value > local_mean)])
  
  quantile_2 = quantile(tribble_scaled_matrix$value[which(tribble_scaled_matrix$cluster == tribble_scaled_matrix$cluster[i] &
                                                            tribble_scaled_matrix$Stade == tribble_scaled_matrix$Stade[i] & 
                                                            tribble_scaled_matrix$value < local_mean)])
  ymin_1 = c(ymin_1, quantile_2[4])
  ymin_2 = c(ymin_2, quantile_2[3])
  ymin_3 = c(ymin_3, quantile_2[2])
  ymin_4 = c(ymin_4, quantile_2[1])
  
  ymax_1 = c(ymax_1, quantile_1[2])
  ymax_2 = c(ymax_2, quantile_1[3])
  ymax_3 = c(ymax_3, quantile_1[4])
  ymax_4 = c(ymax_4, quantile_1[5])
  
}

all_mean = tribble_scaled_matrix %>% group_by(cluster,Stade) %>% summarise(mean = median(value))
all_quantile = tribble_scaled_matrix %>% group_by(cluster,Stade) %>% summarise(q10 = quantile(value, probs = .10),
                                                                                 q20 = quantile(value, probs = .20),
                                                                                 q30 = quantile(value, probs = .30),
                                                                                 q40 = quantile(value, probs = .40),
                                                                                 q50 = quantile(value, probs = .50),
                                                                                 q60 = quantile(value, probs = .60),
                                                                                 q70 = quantile(value, probs = .70),
                                                                               q80 = quantile(value, probs = .80),
                                                                               q90 = quantile(value, probs = .90),
                                                                                 )
for (i in 1:nrow(tribble_scaled_matrix)) {
  
  ymin_1 = c(ymin_1, all_quantile$q40[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                               all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymin_2 = c(ymin_2, all_quantile$q30[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                               all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymin_3 = c(ymin_3, all_quantile$q20[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                               all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymin_4 = c(ymin_4, all_quantile$q10[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                               all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  
  ymax_1 = c(ymax_1, all_quantile$q60[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                              all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymax_2 = c(ymax_2, all_quantile$q70[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                              all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymax_3 = c(ymax_3, all_quantile$q80[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                              all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  ymax_4 = c(ymax_4, all_quantile$q90[which(all_quantile$cluster == tribble_scaled_matrix$cluster[i] &
                                              all_quantile$Stade == tribble_scaled_matrix$Stade[i])])
  
}

tribble_scaled_matrix_test$ymin_1 = ymin_1
tribble_scaled_matrix_test$ymin_2 = ymin_2
tribble_scaled_matrix_test$ymin_3 = ymin_3
tribble_scaled_matrix_test$ymin_4 = ymin_4

tribble_scaled_matrix_test$ymax_1 = ymax_1
tribble_scaled_matrix_test$ymax_2 = ymax_2
tribble_scaled_matrix_test$ymax_3 = ymax_3
tribble_scaled_matrix_test$ymax_4 = ymax_4

RColorBrewer::display.brewer.all()

picked_colors = RColorBrewer::brewer.pal(n = 5,name = "Blues")

final_plot = tribble_scaled_matrix_test %>% ggplot(aes(Stade, value)) +
  geom_line(aes(group = geneID), alpha = 0.3, colour = "lightgrey") +
  geom_line(stat = "summary", fun = "median", colour = "Darkblue", size = 1.5,
            aes(group = 1))+
  geom_ribbon(aes(ymin = ymin_1, ymax = ymax_1, group = 1),
              fill = picked_colors[5], size = 1, linetype = 1, alpha = 0.35)+
  geom_ribbon(aes(ymin = ymin_2, ymax = ymin_1, group = 1),
              fill = picked_colors[4], size = 1, linetype = 1, alpha = 0.35)+
  geom_ribbon(aes(ymin = ymin_3, ymax = ymin_2, group = 1),
              fill = picked_colors[3], size = 1, linetype = 1, alpha = 0.35)+
  geom_ribbon(aes(ymin = ymax_1, ymax = ymax_2, group = 1),
              fill = picked_colors[4], size = 1, linetype = 1, alpha = 0.35)+
  geom_ribbon(aes(ymin = ymax_2, ymax = ymax_3, group = 1),
              fill = picked_colors[3], size = 1, linetype = 1, alpha = 0.35)+
  facet_wrap(facets = vars(cluster), scales = "free") + theme_classic()

save_plots(final_plot,
           width = 14, height= 11,
           graph_name = "Cluster_profile_blue_color")