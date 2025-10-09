# Add annotation to gene clustering output
# Library ####
library(stringr)
library(plyr)

# Parameters ####

# Generated from add_blast_annotation.R
mouse.path = "Mus_musculus_table.tsv"

# Functions ####

add_info = function(current_table, blast_result, specie_name){
  info = c()
  for (i in current_table$`Gene ID`) {
    i = paste(i, "t1", sep = "")
    info = c(info,as.character(blast_result$Name[which(blast_result$TranscriptID == i)]))
  }
  current_table[specie_name] = info
  return(current_table)
}

removequotes = function(list_to_remove){
  return(str_replace_all(list_to_remove, '\"', ""))
}

# Main ####

# Results from expression cluster
file_path = "clustering/hclust_k12"

clusters = 1:12

mouse = read.table(file = mouse.path, sep = "\t", header = T, quote = "")


# Add annotation and create new table. ----

for ( i in clusters) {
  cluster = ifelse(test = i < 10, yes = paste0("0", as.character(i)),no = as.character(i))
  current_table = read.table(file = 
                               paste0(file_path,"/cluster",cluster,"/cluster_gene_list.tsv"),
                             header = F, sep = "\t", quote = "", stringsAsFactors = F)
  colnames(current_table) = "Gene ID"
  current_table = add_info(current_table, mouse, "Mouse")
  write.table(x = current_table,
              file = paste(file_path,"/cluster",cluster,"/cluster_gene_list_annotated.tsv", sep=""),
              quote = F, sep = "\t", row.names = F, col.names = T)
}


# Create Gene_to_cluster_correspondence.tsv file ----


genes = c()
group = c()
count = 0
liste_directories = list.dirs(file_path)
for ( i in liste_directories) {
  if (i != "clustering/hclust_k12") {
    current_table = read.table(file = 
                                 paste0(i,"/cluster_gene_list.tsv"),
                               header = F, sep = "\t", quote = "", stringsAsFactors = F)
    count = count + 1
    genes = c(genes, current_table$V1)
    group = c(group,rep(count, times = nrow(current_table)))
  }
}

df = data.frame(GeneID = genes, Cluster = group)

write.table(x = df, file = "clustering/hclust_k12/cluster_to_gene_correspondance.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)