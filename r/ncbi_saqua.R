library("KEGGREST")
library("keggapi")
library("parallel")
library("data.table")
library("dplyr")

## EMBL ID

# Importar os dados
df_embl_ko <-  read.csv("C:\\Users\\eulle\\Documents\\genkbank_ids\\embl\\uniprot_ko\\df_ko_embl.csv",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
# Criar um chunk dos dados pq são muitas linhas
emb_ko_list <- split(unique(df_embl_ko$KO), ceiling(seq_along(unique(df_embl_ko$KO))/100))

# mapear as vias
list_ko_embl <- list()

for(ii in 1:length(emb_ko_list)){
  kegg_ko_embl <- keggLink("pathway",emb_ko_list[[ii]])
  list_ko_embl <- append(list_ko_embl,kegg_ko_embl)
}

ko_ids_embl <- unique(names(list_ko_embl))
ko_ids_chunk <- split(ko_ids_embl, ceiling(seq_along(ko_ids_embl)/10))

pathway_name_embl <- mapply(KEGGREST::keggGet,ko_ids_embl)

# 27 a 39 é criar uma lista com as vias mapeadas e separá-las com um ;
list_pathway_embl<- list()
sublist_pathway_embl <- list()

for (ii in 1:length(pathway_name)){
  
  sublist_pathway_embl <- c(sublist_pathway_embl,pathway_name[[ii]]$ENTRY)
  sublist_pathway_embl <- append(sublist_pathway_embl,paste(pathway_name[[ii]]$PATHWAY,collapse = ";"))
  list_pathway_embl <-append(list_pathway_embl,list(sublist_pathway_embl))
  sublist_pathway_embl <- list()
  
  
}

# Criar um dataframe com os KOs , uniprot_ids e vias
df_pathnames_embl <- data.frame(Reduce(rbind,list_pathway_embl),stringsAsFactors = FALSE)
rownames(df_pathnames_embl) <- NULL
colnames(df_pathnames_embl) <- c("KO","Pathway")

df_pathway_names_embl <- merge(df_embl_ko, df_pathnames_embl, by="KO", all.x = TRUE)
# Transformar em dataframe para poder salvar em csv

df_embl <- data.frame(lapply(df_pathway_names_embl, as.character), stringsAsFactors=FALSE)
# Tratar os valores duplicados que possuem mais de um ko para uma única proteína
duplicated_embl <- df_embl[duplicated(df_embl$Uniprot_ID) | duplicated(df_embl$Uniprot_ID,fromLast = TRUE),]
df_embl <- df_embl[!(duplicated(df_embl$Uniprot_ID) | duplicated(df_embl$Uniprot_ID,fromLast = TRUE)),]

dat2 <- duplicated_embl %>% group_by(Uniprot_ID) %>% summarise(KO=paste(KO, collapse=","))
dat <- duplicated_embl[!duplicated(duplicated_embl$Pathway) & duplicated_embl$Pathway!="NA",]

duplicated_embl <- merge(dat2,dat,by="Uniprot_ID",all.x = TRUE)
duplicated_embl$KO.y<- NULL
duplicated_embl <-  setnames(duplicated_embl,"KO.x","KO")

df_embl <- rbind(df_embl, duplicated_embl)
rownames(df_embl)<- NULL

#write.csv(df_rev,"pathways_saquaxbiosurfdb_rev",row.names = FALSE,col.names = TRUE)


# REFSEQ ID (Essas linhas abaixo são uma cópia do script acima, exceto pela mudança na linha 111)

df_refseq_ko <-  read.csv("C:\\Users\\eulle\\Documents\\genkbank_ids\\refseq\\uniprot_ko_refseq\\df_ko_refseq.csv",header=TRUE,sep = "\t",stringsAsFactors = FALSE)

refseq_ko_list <- split(unique(df_refseq_ko$KO), ceiling(seq_along(unique(df_refseq_ko$KO))/100))

list_ko_refseq <- list()

for(ii in 1:length(refseq_ko_list)){
  kegg_ko_refseq <- keggLink("pathway",refseq_ko_list[[ii]])
  list_ko_refseq <- append(list_ko_refseq,kegg_ko_refseq)
}

ko_ids_refseq <- unique(names(list_ko_refseq))
ko_ids_chunk_refseq <- split(ko_ids_refseq, ceiling(seq_along(ko_ids_refseq)/10))

pathway_name_refseq <- mapply(KEGGREST::keggGet,ko_ids_refseq)


list_pathway_refseq<- list()
sublist_pathway_refseq <- list()

for (ii in 1:length(pathway_name_refseq)){
  
  sublist_pathway_refseq <- c(sublist_pathway_refseq,pathway_name_refseq[[ii]]$ENTRY)
  sublist_pathway_refseq <- append(sublist_pathway_refseq,paste(pathway_name_refseq[[ii]]$PATHWAY,collapse = ";"))
  list_pathway_refseq <-append(list_pathway_refseq,list(sublist_pathway_refseq))
  sublist_pathway_refseq <- list()
  
  
}

df_pathnames_refseq <- data.frame(Reduce(rbind,list_pathway_refseq),stringsAsFactors = FALSE)
rownames(df_pathnames_refseq) <- NULL
colnames(df_pathnames_refseq) <- c("KO","Pathway")

df_pathway_names_refseq <- merge(df_refseq_ko, df_pathnames_refseq, by="KO", all.x = TRUE)
# Transformar em dataframe para poder salvar em csv
df_refseq <- data.frame(lapply(df_pathway_names_refseq, as.character), stringsAsFactors=FALSE)

duplicated_refseq <- df_refseq[duplicated(df_refseq$Uniprot_ID) | duplicated(df_refseq$Uniprot_ID,fromLast = TRUE),]
df_refseq <- df_refseq[!(duplicated(df_refseq$Uniprot_ID) | duplicated(df_refseq$Uniprot_ID,fromLast = TRUE)),]

dat2_refseq <- duplicated_refseq %>% group_by(Uniprot_ID) %>% summarise(KO=paste(KO, collapse=","))
dat_refseq <- duplicated_refseq[duplicated(duplicated_refseq$Uniprot_ID),]

duplicated_refseq <- merge(dat2_refseq,dat_refseq,by="Uniprot_ID",all.x = TRUE)
duplicated_refseq$KO.y<- NULL
duplicated_refseq <-  setnames(duplicated_refseq,"KO.x","KO")

df_refseq <- rbind(df_refseq, duplicated_refseq)
rownames(df_refseq)<- NULL

# Concatenar os embl_ids e refseq_ids

df_genkbank <- rbind(df_refseq,df_embl)
df_genkbank <- df_genkbank[!duplicated(df_genkbank$Uniprot_ID),]

write.csv(df_genkbank,"pathways_ncbixsaqua",row.names = FALSE,col.names = TRUE)
