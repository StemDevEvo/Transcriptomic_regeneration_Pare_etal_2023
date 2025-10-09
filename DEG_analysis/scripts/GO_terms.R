# Script to generate GOterms annotation / enrichment from a list of gene.

library(AnnotationHub)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)

### Function


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

match_id = function(id, blast_table){
  return(blast_table$BlastResult[which(blast_table$TranscriptID == paste(id,"t1",sep=""))])
}


####



#########################

# OrgDb : Mouse.

hub <- AnnotationHub()

q <- query(hub, "org.Mm.eg.db")
id <- q$ah_id[length(q)]
OrgDb = hub[[id]]

##########################

# Loading TMM matrices to determine universe (all genes initially compared)

matrix_NA_1TMM = read_csv(file = "Table/NA/NA_1TMM_threshold_matrix.csv")
matrix_0dpa_1TMM = read_csv(file = "Table/0dpa/0dpa_1TMM_threshold_matrix.csv")
matrix_1dpa_1TMM = read_csv(file = "Table/1dpa/1dpa_1TMM_threshold_matrix.csv")
matrix_2dpa_1TMM = read_csv(file = "Table/2dpa/2dpa_1TMM_threshold_matrix.csv")
matrix_3dpa_1TMM = read_csv(file = "Table/3dpa/3dpa_1TMM_threshold_matrix.csv")
matrix_5dpa_1TMM = read_csv(file = "Table/5dpa/5dpa_1TMM_threshold_matrix.csv")


universe = unique(c(matrix_NA_1TMM$...1,
                    matrix_0dpa_1TMM$...1,
                    matrix_1dpa_1TMM$...1,
                    matrix_2dpa_1TMM$...1,
                    matrix_3dpa_1TMM$...1,
                    matrix_5dpa_1TMM$...1))

###########################

# Load mouse annotation :

ann_table = read.table(file = "Mus_musculus_table.tsv", 
                       sep = "\t", header = T)


##########################

# Universe ----

universe_specie = lapply(X = universe, FUN = match_id, blast_table = ann_table)

df = data.frame(universe_transcript = universe,
                universe_specie = unlist(universe_specie))

universe_gene = df$universe_specie[df$universe_specie != "."]
universe_gene = unlist(lapply(X = universe_gene, FUN = return_uniprot_id))

# List to compare ----

liste_directories = list.dirs(path = "DEG")

for (i in liste_directories){
  if (i != "DEG") {
    i = strsplit(i, split = "/", fixed = T)[[1]][length(strsplit(i, split = "/", fixed = T)[[1]])]
    DGE_result = read_tsv(file = paste0("DEG/",i,"/DGE_result_FC_and_FDR_filter_fully_annotated.tsv"))
    
    
    DGE_gene = unlist(lapply(DGE_result$Mouse, prepare_DGE_gene))
    
    DGE_gene = unlist(lapply(X = DGE_gene, FUN = return_uniprot_id))
    
    
    ## GO terms enrichment (one condition only)
    
    file_path = paste("GOterms/",i,"/",sep ="")
    if(!dir.exists(file_path)){
      dir.create(file_path)
    }
    ont = c("BP","MF", "CC")
    for (y in ont) {
      ego <- enrichGO(gene = DGE_gene,
                      universe = universe_gene,
                      OrgDb = OrgDb,
                      keyType = "UNIPROT",
                      ont = y,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = T)
      
      write.table(x = ego, file = paste(file_path,"GO_terms_enrichment_",y,".tsv", sep = ""), sep = "\t",quote = F)
      
      
      if(length(which(ego@result$p.adjust < 0.05)) > 5){
        
        pdf(file = paste(file_path,"GO_terms_plot_",y,".pdf",sep=""))
        print(goplot(ego, showCategory = 10))
        dev.off()
        
        pdf(file = paste(file_path,"barplot_top20_",y,".pdf",sep=""))
        print(barplot(ego, showCategory = 20))
        dev.off()
        
        pdf(file = paste(file_path,"dotplot_top20_",y,".pdf", sep =""))
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        
        egox = setReadable(ego, OrgDb, "UNIPROT")
        
        egox2 = pairwise_termsim(egox)
        pdf(file = paste(file_path,"Tree_plot_clustering_",y,".pdf",sep=""), width = 17, height = 8)
        print(treeplot(egox2, hclust_method = "average"))
        dev.off()
        
        ego_enr = pairwise_termsim(ego)
        
        pdf(file = paste(file_path,"Network_plot",y,".pdf",sep=""), width = 12, height = 12)
        print(emapplot(ego_enr))
        dev.off()
      }
    }
  }
}


################# Kegg pathways
library(pathview)

get_logFC = function(DGE_result){
  if(DGE_result$Mouse == "."){
    return()
  }else{
    return(DGE_result$X.logFC.)
  }
}

return_gene_symbol <- function(uniprot){
  return(strsplit(x = str_extract(uniprot, "GN=[^ ]*"), split = "=", fixed = T)[[1]][2])
}

for (i in liste_directories){
  cond = strsplit(i, split = "/", fixed = T)[[1]][length(strsplit(i, split = "/", fixed = T)[[1]])]
  if (i != "DEG") {
    if(!dir.exists(paste0("Kegg/",cond))){
      dir.create(paste0("Kegg/",cond))
    }
    file_path = paste("Kegg/",cond,"/",sep ="")
    DGE_result = read_tsv(file = paste0(i,"/DGE_result_FC_and_FDR_filter_fully_annotated.tsv"))
    DGE_gene = unlist(lapply(DGE_result$Mouse, prepare_DGE_gene))
    DGE_gene = unlist(lapply(X = DGE_gene, FUN = return_uniprot_id))
    
    #DGE_gene = unlist(lapply(X = DGE_result$Mouse, FUN = return_gene_symbol))
    
    logFC = ifelse(DGE_result$Mouse ==".",NA,DGE_result$X.logFC.)
    logFC = logFC[!is.na(logFC)]
    names(logFC) = DGE_gene
    kk <- enrichKEGG(gene = DGE_gene,
                     universe = universe_gene,
                     keyType = "uniprot",
                     organism = "mmu",
                     pvalueCutoff = 0.05)
    write.table(x = kk, file = paste0(file_path,"Kegg_pathways_enrichment.tsv"), sep = "\t",quote = F, col.names = NA)
  
    
    ENTREZ_name_logfc = id2eg(ids = names(logFC), category = "UNIPROT", org = "Mm")
    
    names(logFC) = ENTREZ_name_logfc[,2]
    
    significant_ids = kk@result$ID[1:length(which(kk@result$p.adjust < 0.05))]
    
    for (y in 1:5) {
      current_directory = getwd()
      setwd(file_path)
      current_pathway_id = significant_ids[y]
      if(!is.na(current_pathway_id)){
        pathview(gene.data = logFC,
               pathway.id = current_pathway_id,
               gene.idtype = "ENTREZ",
               species = "mmu",
               limit = list(gene=max(abs(logFC)), cpd=1))
        setwd(current_directory)
      }
      setwd(current_directory)
    }
    setwd(current_directory)
  }
}
