library("KEGGREST")
#library("keggapi")
library("tidyverse")


df_ko = read.csv("C:\\Users\\eulle\\Documents\\uniprot-ko_id\\CoSQG/df_ko_ids_CoSQG.csv",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
ko_rev <- KEGGREST::keggLink("pathway",unique(df_ko$KO))
ko_ids_rev <- unique(names(ko_rev))

# Saber o que significa ceiling e set_along
# Fiz o chunk porque o keggget só aceita uma lista de no max 10 itens
ko_ids_rev <- split(ko_ids_rev, ceiling(seq_along(ko_ids_rev)/10))


pathway_name_rev <- mapply(KEGGREST::keggGet,ko_ids_rev)
pathway_name_rev <- unlist(pathway_name_rev,recursive = FALSE)
# extraindo as vias metabólicas capturadas no KEGGREST
list_pathway_rev <- list()
sublist_pathway_rev <- list()
dict_pathway <- list()
dict_aux <- list()

for (ii in 1:length(pathway_name_rev)){
  
  
  sublist_pathway_rev <- c(sublist_pathway_rev,pathway_name_rev[[ii]]$ENTRY)
  dict_aux <- as.list(pathway_name_rev[[ii]]$PATHWAY)
  dict_aux <- Biobase::reverseSplit(dict_aux)
  dict_pathway <- append(dict_aux,dict_pathway)
  dict_pathway <- dict_pathway[!duplicated(dict_pathway)]
  sublist_pathway_rev <- append(sublist_pathway_rev,paste(pathway_name_rev[[ii]]$PATHWAY,collapse = ";"))
  list_pathway_rev <-append(list_pathway_rev,list(sublist_pathway_rev))
  sublist_pathway_rev <- list()
  
  
  
}



# Juntar list_pathway e saqua_biosurfdb

df_pathnames_rev <- data.frame(Reduce(rbind,list_pathway_rev),stringsAsFactors = FALSE)
rownames(df_pathnames_rev) <- NULL
colnames(df_pathnames_rev) <- c("KO","Pathway")

df_pathway_names_rev <- merge(df_ko, df_pathnames_rev, by="KO", all.x = TRUE)

df_rev <- data.frame(lapply(df_pathway_names_rev, as.character), stringsAsFactors=FALSE)
df_rev <- df_rev[!duplicated(df_rev[,"Uniprot_ID"]),]

write.csv(df_rev,"pathways_biosufdbxROCHA2",row.names = FALSE,col.names = TRUE)

# Fazer uma lista de KO para cada rota.

# Verificar se há uniprot duplicado
df_kos_count <- data.frame(matrix(ncol = 6, nrow = 0),stringsAsFactors = FALSE)
x <- c("Name","cod_pathway", "KO", "Count","KO_found","KO_per_total")
colnames(df_kos_count) <- x
unique_ko <- unique(df_ko$KO) # Quantos foram os KOS unicos foram mapeados na análise
for(index in 1:length(dict_pathway)){
  ko_p_rota <- keggLink("ko", unname(dict_pathway[index]))
  kos <- gsub("^.*\\:(.*)", "\\1", unname(ko_p_rota))
  # Intersecção dos dos ko
  count_ko <- length(kos)
  kos_concat <- paste(kos,collapse = ";")
  df_kos_count <- df_kos_count %>% add_row(Name = names(dict_pathway[index]),cod_pathway=unname(dict_pathway[index]), KO = kos_concat, Count = count_ko,KO_found = paste(intersect(unique_ko,kos),collapse = ";"), KO_per_total = (length(intersect(unique_ko,kos))/count_ko)*100)
  df_kos_count[,"KO_per_total"] <- round(df_kos_count[,"KO_per_total"],2)
}

df_kos_count<- data.frame(lapply(df_kos_count[,c("Name","cod_pathway","KO","KO_found")], as.character),Count=as.integer(df_kos_count$Count),KO_per_total=as.numeric(df_kos_count$KO_per_total), stringsAsFactors=FALSE)
write.table(df_kos_count,"pathways_countsxCoSQG.csv",row.names = FALSE,col.names = TRUE)



