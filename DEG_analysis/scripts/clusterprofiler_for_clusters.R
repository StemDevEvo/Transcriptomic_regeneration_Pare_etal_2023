# Apply cluterProfiler to gene expression clustering

library(clusterProfiler)
library(AnnotationHub)
library(enrichplot)

## Parameters :

cutoff_char = "_pval_cutoff"


### Function :

match_id = function(id, blast_table){
  return(blast_table$BlastResult[which(blast_table$TranscriptID == paste(id,"t1",sep=""))])
}

prepare_DGE_gene = function(blast_res){
  if(blast_res == "."){
    return()
  }else{
    return(unlist(strsplit(blast_res, " "))[1])
  }
}

return_uniprot_id = function(full_name){
  splitted_id = unlist(strsplit(x=full_name, split = "|", fixed = T))
  return(splitted_id[2])
}

#########################

# OrgDb : Mouse.

hub <- AnnotationHub()

q <- query(hub, "org.Mm.eg.db")
id <- q$ah_id[length(q)]
OrgDb = hub[[id]]

##########################

# Loading TMM matrices to determine universe (all genes initially compared)

NA_1tmm = read.table("Table/NA/NA_1TMM_threshold_matrix.csv",
                     header = T, row.names = 1, com='', sep = ",")

dpa0_1tmm = read.table("Table/0dpa/0dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")

dpa1_1tmm = read.table("Table/1dpa/1dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")

dpa2_1tmm = read.table("Table/2dpa/2dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")

dpa3_1tmm = read.table("Table/3dpa/3dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")

dpa5_1tmm = read.table("Table/5dpa/5dpa_1TMM_threshold_matrix.csv",
                       header = T, row.names = 1, com='', sep = ",")



# Correspond to all the genes tested for DGE

universe = unique(c(rownames(NA_1tmm),rownames(dpa0_1tmm),
                    rownames(dpa1_1tmm),rownames(dpa2_1tmm),
                    rownames(dpa3_1tmm),rownames(dpa5_1tmm)))

###########################

# Load human or mouse annotation : (see add_blast_annotation.R)

ann_table = read.table(file = "Mus_musculus_table.tsv", 
                       sep = "\t", header = T)

##########################

# Determine GO terms enrichment for each cluster gene list.

clusters = 1:12

# Take the universe variable and translate the geneID from Platynereis to Mmus gene ID.
universe_specie = lapply(X = universe, FUN = match_id, blast_table = ann_table)

for (i in clusters){
  name = ifelse(i < 10, paste("0",as.character(i),sep=""),i)
  # List of gene with Platynereis geneID present in a given cluster
  # one geneID per line
  gene_list = read.table(
  file = paste("clustering/hclust_k12/cluster",name,"/cluster_gene_list.tsv",sep = ""),
  sep = "\t", header = F)
  
  # Translate the gene id list to Mmus
  unordered_list_mouse = unlist(lapply(X = gene_list$V1, FUN = match_id, blast_table = ann_table))
  

  # formatting data
  df = data.frame(universe_transcript = universe,
                  universe_specie = unlist(universe_specie))
  DGE_gene = unlist(lapply(unordered_list_mouse, prepare_DGE_gene))
  DGE_gene = unlist(lapply(X = DGE_gene, FUN = return_uniprot_id))
  
  ## GO terms enrichment (one condition only)
  
  universe_gene = df$universe_specie[df$universe_specie != "."]
  universe_gene = unlist(lapply(X = universe_gene, FUN = return_uniprot_id))
  ont = c("BP","MF", "CC")
  for (y in ont) {
    ego <- enrichGO(gene = DGE_gene,
                    universe = universe_gene,
                    OrgDb = OrgDb,
                    keyType = "UNIPROT",
                    ont = y,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = T
    )
    
    
    file_path = paste("clustering/hclust_k12/cluster",name,"/",sep ="")
    
    write.table(x = ego@result, file = paste(file_path,"GO_terms_enrichment_",y,cutoff_char,".tsv", sep = ""),
                sep = "\t",quote = F, col.names = NA)
    
    
    if(length(which(ego@result$p.adjust < 0.05)) > 5){
      ego@result$qvalue = ego@result$p.adjust
      
      pdf(file = paste(file_path,"GO_terms_plot_",y,cutoff_char,".pdf",sep=""))
      print(goplot(ego, showCategory = 10))
      dev.off()
      
      pdf(file = paste(file_path,"barplot_top20_",y,cutoff_char,".pdf",sep=""))
      print(barplot(ego, showCategory = 20))
      dev.off()
      
      pdf(file = paste(file_path,"dotplot_top20_",y,cutoff_char,".pdf", sep =""))
      print(dotplot(ego, showCategory = 20))
      dev.off()
      
      
      egox = setReadable(ego, OrgDb, "UNIPROT")
      
      egox2 = pairwise_termsim(egox)
      pdf(file = paste(file_path,"Tree_plot_clustering_",y,cutoff_char,".pdf",sep=""), width = 17, height = 8)
      print(treeplot(egox2, hclust_method = "average"))
      dev.off()
      
      ego_enr = pairwise_termsim(ego)
      
      pdf(file = paste(file_path,"Network_plot_",y,cutoff_char,".pdf",sep=""), width = 12, height = 12)
      print(emapplot(ego_enr))
      dev.off()
    }
  }
}