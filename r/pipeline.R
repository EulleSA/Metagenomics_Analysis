#!/usr/bin/env Rscript

# conda install -c conda-forge r-devtools

#update.packages(ask = FALSE,
#                checkBuilt = TRUE,repos="http://cran.wustl.edu/")

packages = c("future", "stringr","remotes")

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager",repos="http://cran.wustl.edu/")
    BiocManager::install("KEGGREST",type="source")
    BiocManager::install("Biobase")
}

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x,repos="http://cran.wustl.edu/",type="source")
      library(x, character.only=TRUE)
    }
  }
)

#install.packages("furrr",upgrade=FALSE,repos="http://cran.wustl.edu/")

#remove.packages("KEGGREST")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")


if(!require("tidyverse")){
    devtools::install_github("tidyverse/tidyverse")
}

if(!require("thacklr")){
  
  devtools::install_github("thackl/thacklr",auth_token = "dcfbe5870d9e9e463073938c453105d74a0e5ef7",force=TRUE,upgrade=FALSE)
  
}



library(stringr)
library(future)
library(furrr)
library("tidyverse")
library("thacklr")
library(KEGGREST)
library(purrr)
# colnames_eggnog <- "query_name.seed_eggNOG_ortholog.seed_ortholog_evalue.seed_ortholog_score.best_tax_level.Preferred_name.GOs.EC.KEGG_ko.KEGG_Pathway.KEGG_Module.KEGG_Reaction.KEGG_rclass.BRITE.KEGG_TC"
# colnames_eggnog <- strsplit(colnames_eggnog,split = "\\.")[[1]]
# colnames(df_eggnog) <- colnames_eggnog

#df_eggnog_2 <- read.csv("D:/Downloads/CoSqG.emapper.csv",sep=",",stringsAsFactors = FALSE)
#paste(strsplit(colnames(df_eggnog_2),split = "\\.")[[1]][2:16])

args = commandArgs(trailingOnly = TRUE)

dataframe_final <- function(df_eggnog_csv,file_name){
  
  file_name = strsplit(strsplit(file_name,split = "/")[[1]][3],split="\\.")[[1]][1]
  
  #names(df_eggnog_csv)[names(df_eggnog_csv) == "X.query_name"] <- "query_name"
  #df_eggnog_csv <- subset(df_eggnog_csv,KEGG_ko != "",select = c("query_name","best_tax_level","KEGG_ko","KEGG_Pathway"))
  colnames_eggnog <- "query_name.seed_eggNOG_ortholog.seed_ortholog_evalue.seed_ortholog_score.best_tax_level.Preferred_name.GOs.EC.KEGG_ko.KEGG_Pathway.KEGG_Module.KEGG_Reaction.KEGG_rclass.BRITE.KEGG_TC"
  colnames_eggnog <- strsplit(colnames_eggnog,split = "\\.")[[1]]
  colnames(df_eggnog_csv) <- colnames_eggnog
  
  df_eggnog_csv$KEGG_ko <- gsub("ko:","",df_eggnog_csv$KEGG_ko)
  df_eggnog_csv$KEGG_Pathway <- gsub(",map[0-9]+","",df_eggnog_csv$KEGG_Pathway)
  
  df_eggnog_csv <- subset(df_eggnog_csv,KEGG_ko != "",select = c("query_name","best_tax_level","KEGG_ko","KEGG_Pathway"))
  
  pathways_names_kos <- unique(unlist(str_split(df_eggnog_csv[is.na(df_eggnog_csv$KEGG_Pathway) == FALSE,"KEGG_ko"][[1]],","),use.names = FALSE))
  kos_pathways_chunk <- split(pathways_names_kos, ceiling(seq_along(pathways_names_kos)/10))
  
  # pathways_cod <- unique(unlist(str_split(df_eggnog_csv[df_eggnog_csv$KEGG_Pathway != "","KEGG_Pathway"],","),use.names = FALSE))
  # cod_pathways_chunk <- split(pathways_cod, ceiling(seq_along(pathways_cod)/10))
  # 
  # 
  # codpathways_name <- lapply(cod_pathways_chunk,KEGGREST::keggGet)
 
  numCores <- parallel::detectCores() - 1
  print(parallel::detectCores())
  plan(multisession, workers = numCores)
  
  kospathways_name <- furrr::future_map(kos_pathways_chunk,keggGet)
  kospathways_name <- unlist(kospathways_name,recursive = FALSE)
  
  # extraindo as vias metabólicas capturadas no KEGGREST
  list_pathway_rev <- list()
  sublist_pathway_rev <- list()
  dict_pathway <- list()
  dict_aux <- list()
  
  for (ii in 1:length(kospathways_name)){
    
    
    sublist_pathway_rev <- c(sublist_pathway_rev,kospathways_name[[ii]]$ENTRY)
    dict_aux <- as.list(kospathways_name[[ii]]$PATHWAY)
    dict_aux <- Biobase::reverseSplit(dict_aux)
    dict_pathway <- append(dict_aux,dict_pathway)
    dict_pathway <- dict_pathway[!duplicated(dict_pathway)]
    sublist_pathway_rev <- append(sublist_pathway_rev,paste(kospathways_name[[ii]]$PATHWAY,collapse = ";"))
    list_pathway_rev <-append(list_pathway_rev,list(sublist_pathway_rev))
    sublist_pathway_rev <- list()
    
  }
  
  df_pathnames_rev <- data.frame(Reduce(rbind,list_pathway_rev),stringsAsFactors = FALSE)
  
  rownames(df_pathnames_rev) <- NULL
  colnames(df_pathnames_rev) <- c("KO","Pathway")
  
  df_kos_count <- data.frame(matrix(ncol = 6, nrow = 0),stringsAsFactors = FALSE)
  x <- c("Name","cod_pathway", "KO","KO_found", "Count","KO_per_total")
  colnames(df_kos_count) <- x
  unique_ko <- unique(pathways_names_kos) # Quantos foram os KOS unicos foram mapeados na análise
  for(index in 1:length(dict_pathway)){
    ko_p_rota <- KEGGREST::keggLink("ko", unname(dict_pathway[index]))
    kos <- gsub("^.*\\:(.*)", "\\1", unname(ko_p_rota))
    # Interseccao dos dos ko
    count_ko <- length(kos)
    kos_concat <- paste(kos,collapse = ";")
    df_kos_count_row <- dplyr::bind_rows(Name = gsub("/",":",names(dict_pathway[index])),cod_pathway=unname(dict_pathway[index])[[1]], KO = kos_concat,KO_found = paste(intersect(unique_ko,kos),collapse = ";"),Count = length(intersect(unique_ko,kos)), KO_per_total = (length(intersect(unique_ko,kos))/count_ko)*100)
    df_kos_count <- rbind(df_kos_count,df_kos_count_row)
    df_kos_count[,"KO_per_total"] <- round(df_kos_count[,"KO_per_total"],2)
  }
  write.table(df_kos_count,paste0("results/",file_name,"_count_pathways"),sep="\t",row.names = FALSE,col.names = TRUE)
  write.table(na.omit(df_eggnog_csv[df_eggnog_csv$KEGG_Pathway!="NA","query_name"]),paste0("results/",file_name,"_reads_ids.txt"),row.names = FALSE,col.names = FALSE,quote = FALSE)
  #return(df_kos_count)
}

df_eggnog <- read_eggnog(args[1])

dataframe_final(df_eggnog,args[1])

#df_kos_count_final <- dataframe_final(df_eggnog)

#write.table(df_kos_count_final,"D:/Downloads/pathways_count_CoSqG",sep="\t",row.names = FALSE,col.names = TRUE)

