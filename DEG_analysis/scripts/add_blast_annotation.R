library(stringr)
library(plyr)

# Results from blasts (Transcriptome to proteome of selected species)
mouse.path = "Mus_musculus_table.tsv"
human.path = "Hsap_table.tsv"
aplysia.path = "Aplysia_californica_table.tsv"
drosophila.path = "R_Drosophila_table.tsv"
saccoglossus.path = "R_Saccoglossus_kowalevskii_table.tsv"
zebrafish.path = "R_Danio_rerio.tsv"

add_info = function(current_table, blast_result, specie_name){
  info = c()
  for (i in current_table$X..) {
    i = paste(i, "t1", sep = "")
    info = c(info,as.character(blast_result$Name[which(blast_result$TranscriptID == i)]))
  }
  current_table[specie_name] = info
  return(current_table)
}

removequotes = function(list_to_remove){
  return(str_replace_all(list_to_remove, '\"', ""))
}


liste_directories = list.dirs("analysis/DEG2")

mouse = read.table(file = mouse.path, sep = "\t", header = T, quote = "")
human = read.table(file = human.path, sep = "\t", header = T, quote = "")
aplysia = read.table(file = aplysia.path, sep = "\t", header = T, quote = "")
drosophila = read.table(file = drosophila.path, sep = "\t", header = T, quote = "")
saccoglossus = read.table(file = saccoglossus.path, sep = "\t", header = T, quote = "")
zebrafish = read.table(file = zebrafish.path, sep = "\t", header = T, quote = "")


for ( i in liste_directories) {
  if (i != "DEG") {
    current_table = read.table(file = 
                                 paste0(i,"/DGE_result_",
                                        strsplit(i, split = "/")[[1]][3],
                                        "_FC_and_FDR_filter_trinotate.tsv"),
                               header = T, sep = "\t", quote = "", stringsAsFactors = F)
    current_table$X.. = removequotes(current_table$X..)
    rownames(current_table) = current_table$X..
    current_table = add_info(current_table, mouse, "Mouse")
    current_table = add_info(current_table, human, "Human")
    current_table = add_info(current_table, drosophila, "Drosophila")
    current_table = add_info(current_table, aplysia, "Aplysia californica")
    current_table = add_info(current_table, saccoglossus, "Saccoglossus kowalevskii")
    current_table = add_info(current_table, zebrafish, "Danio rerio")
    write.table(x = current_table,
                file = paste(i,"/DGE_result_FC_and_FDR_filter_fully_annotated.tsv", sep=""),
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



for (i in liste_directories) {
  if (i != ".") {
    table_5species = read.table(file = paste(i,"/DGE_result_FC_and_FDR_filter_annotated_5species.tsv", sep=""),
                                header = T, sep = "\t", quote = "", stringsAsFactors = F)
    table_annotated = read.table(file = paste(i,"/copie.tsv", sep=""),
                                 header = T, sep = "\t", quote = "", stringsAsFactors = F, strip.white = T)
    table_annotated["Mouse"] = table_5species$Mouse
    table_annotated["Human"] = table_5species$Human
    table_annotated["Drosophila"] = table_5species$Drosophila
    table_annotated["Aplysia clifornica"] = table_5species$Aplysia.clifornica
    table_annotated["Saccoglossus kowalevskii"] = table_5species$Saccoglossus.kowalevskii
    write.table(x = table_annotated,
                file = paste(i,"/DGE_result_FC_and_FDR_filter_annotated_fully_annotated.tsv", sep=""),
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



