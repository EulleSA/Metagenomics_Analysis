#!/usr/bin/env Rscript



library(KEGGREST)
library(KEGGgraph)
library(igraph)
library(tidyr)
library(dplyr)

library(ggraph)
library(scatterpie)
library(easylayout)
library(stringr)
library(ggiraph)

args = commandArgs(trailingOnly=TRUE)

#df_pathway <- df_pathway %>% add_row(node1 = "K00529", node2 = "K00529")
pathway2dataframe <- function(cod_pathway){
  kgml <- suppressMessages(KEGGREST::keggGet(cod_pathway, "kgml"))
  mapkpathway <- KEGGgraph::parseKGML(kgml)
  mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway,FALSE)
  aux <- names(mapkG@edgeData@data)
  aux <- as.data.frame(aux, stringsAsFactors = FALSE)
  aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
  colnames(aux)[1] <- "node1"
  aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)
  
  nodes <- nodes(mapkG)[grepl("^ko",nodes(mapkG))]
  nodes_noEdge <- nodes[!(nodes %in% aux$node1 | nodes %in% aux$node2)]
  
  if(length(nodes_noEdge)!=0)
  {
    nodes_noEdge <- lapply(nodes_noEdge,rbind)
    nodes_noEdge <- cbind(nodes_noEdge,nodes_noEdge)
    colnames(nodes_noEdge) <- c("node1","node2")
  }
  
  rm(mapkpathway, kgml,nodes)
  aux <- rbind(aux,nodes_noEdge)
  return(aux)

}

uncollapse <- function(df_pathway){
  colnames(df_pathway) <- c("key", "value")
  df_pathway <- df_pathway %>%
    mutate(key = strsplit(key, ";")) %>%
    unnest(key)
  return(df_pathway)
}
# caso de uso


df_pathway <- pathway2dataframe(args)

idx <- grepl("^path", df_pathway$node1) | grepl("^path", df_pathway$node2)
df_pathway <- df_pathway[!idx,]

df_pathway$node1 <- gsub("^ko:", "", df_pathway$node1)
df_pathway$node2 <- gsub("^ko:", "", df_pathway$node2)

df_pathway <- unique(uncollapse(df_pathway))
df_pathway <- unique(uncollapse(df_pathway[, c(2,1)]))[, c(2,1)]
colnames(df_pathway) <- c("node1", "node2")


graph_kegg <- graph_from_data_frame(df_pathway, TRUE)


# Construction of the graphs pathways

# Criar um dataframe
complete_dataframe <- function(rowname_taxonomy,cod_pathway){
  
  #df_kos_count <- 
  #rowname_taxonomy <- df_taxonomybesthits[unique(unlist(str_split(df_taxonomybesthits$KEGG_ko,","))) %in% str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]] ,c("KEGG_ko","best_tax_level")]
  rowname_taxonomy <- df_taxonomybesthits[unlist(lapply(str_split(df_kos_count_final[df_kos_count_final$cod_pathway == cod_pathway,"KO_found"],";")[[1]], function(x) grep(x, df_taxonomybesthits$KEGG_ko, fixed = TRUE))),c("KEGG_ko","best_tax_level")]
  
  
  name_org <- unique(rowname_taxonomy$best_tax_level)
  colname_kos <- V(graph_kegg)$name
  kos_binary <- lapply(name_org,function(name_org){as.numeric(colname_kos %in% unique(unlist(str_split(rowname_taxonomy[rowname_taxonomy$best_tax_level==name_org,"KEGG_ko"],","))))})
  
  if(length(kos_binary) == 0){
    df_categories <- data.frame(matrix(ncol = length(colname_kos), nrow = 0))
    colnames(df_categories) <- colname_kos
    
    `%notin%`<- Negate(`%in%`)
    unclassified_kos <- str_split(df_kos_count_final[df_kos_count_final$cod_pathway == cod_pathway,"KO_found"],";")[[1]]
    unclassified_kos <- unclassified_kos[unclassified_kos %notin% unique(unlist(str_split(rowname_taxonomy$KEGG_ko,",")))]
    
    df_categories <- rbind(df_categories,ifelse(colname_kos %in% unclassified_kos,1,0))
    rownames(df_categories)[nrow(df_categories)] <- "Unclassified"
    
    df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>=1,0,1))
    rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  }
  else{
    df_categories <- data.frame(matrix(unlist(kos_binary), nrow=length(kos_binary), byrow=T),stringsAsFactors=FALSE)
    
    `%notin%`<- Negate(`%in%`)
    unclassified_kos <- str_split(df_kos_count_final[df_kos_count_final$cod_pathway == cod_pathway,"KO_found"],";")[[1]]
    unclassified_kos <- unclassified_kos[unclassified_kos %notin% unique(unlist(str_split(rowname_taxonomy$KEGG_ko,",")))]
    
    rownames(df_categories) <- name_org
    colnames(df_categories) <- V(graph_kegg)$name
    
    df_categories <- rbind(df_categories,ifelse(colnames(df_categories) %in% unclassified_kos,1,0))
    rownames(df_categories)[nrow(df_categories)] <- "Unclassified"
    
    df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>=1,0,1))
    rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  }
  
  
  return(df_categories)
}

df_taxonomybesthits <- read.csv("D:/Downloads/taxmy_CoSqG",sep = "\t",stringsAsFactors = FALSE)
df_kos_count_final <- read.csv("D:/Downloads/pathways_count_CoSqG",sep='\t', stringsAsFactors = FALSE)

# Gerando layout

layout <- easylayout::vivagraph(graph_kegg)


df_organisms <- complete_dataframe(df_taxonomybesthits,args)

# Inserindo dados nos vertices
# Categoria é um valor binário onde 1 - COR e 0 - NÃO COR
# x e y são as coordenadas
# Ver esses negócios dA LISTA
categories <- list(name=V(graph_kegg)$name)

categories <- append(categories,lapply(lapply(setNames(split(df_organisms,seq(nrow(df_organisms))),rownames(df_organisms)),unname),unlist))
categories[["x"]] <- layout[,1]
categories[["y"]] <- layout[,2]
vertex_attr(graph_kegg) <- categories

# Inserindo o tamanho dos nós
# PERGUNTAR COMO COLOCA O MESMO TAMANHO PARA TODOS 
#V(graph_kegg)$size <- V(graph_kegg)$categoria 

svg(filename=paste(args,".svg"),wi)
# Plotando
ggraph(graph_kegg, x = x, y = y) +
  
  # Arestas primeiro
  geom_edge_link(
    aes(end_cap = circle(1.5 + 0.5, 'mm'))
    ,color = "#808080"
    ,arrow = arrow(length = unit(1.5, 'mm'), type = "closed")
    ,alpha = 0.5
  ) +
  
  # Pie charts
  geom_scatterpie(
    data    = igraph::as_data_frame(graph_kegg, "vertices")
    ,mapping = aes(x, y, group = name, r = 7)
    ,cols    = rownames(df_organisms)
    ,color   = "black"
  ) +
  
  # Aspect ratio
  coord_fixed() +
  # Node label
  geom_node_text(mapping = aes(x, y), label = gsub("^ec:", "", V(graph_kegg)$name),
                 nudge_y = 15) +
  theme_void()

dev.off()

