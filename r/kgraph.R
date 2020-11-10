library(KEGGREST)
library(KEGGgraph)
library(igraph)
library(tidyr)
library(dplyr)

# pathway2dataframe ####

#' Get the edges from a given KEGG pathway
#'
#' Given a KEGG pathway ID, this function returns a data.frame ready to create
#' an igraph object.
#'
#' @param pathway A KEGG pathway ID.
#' 
#' @return This function returns a data.frame containing the edges from a
#' KEGG pathway.
#'
#' @examples
#' \dontrun{
#' df <- pathway2dataframe("ko00010")
#' }
#'
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGgraph parseKGML
#' @importFrom KEGGgraph KEGGpathway2Graph
#' @importFrom KEGGREST keggLink
#'
#' @author
#' Diego Morais

#df_pathway <- df_pathway %>% add_row(node1 = "K00529", node2 = "K00529")
pathway2dataframe <- function(pathway){
  kgml <- suppressMessages(KEGGREST::keggGet(pathway, "kgml"))
  mapkpathway <- KEGGgraph::parseKGML(kgml)
  mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway,FALSE)
  nodes <- nodes(mapkG)[grepl("^ko",nodes(mapkG))]
  nodes_noEdge <- nodes[!(nodes %in% df_pathway$node1 | nodes %in% df_pathway$node2)]
  nodes_noEdge <- lapply(nodes_noEdge,rbind)
  nodes_noEdge <- cbind(nodes_noEdge,nodes_noEdge)
  colnames(nodes_noEdge) <- c("node1","node2")
  rm(mapkpathway, kgml,nodes)
  aux <- names(mapkG@edgeData@data)
  aux <- as.data.frame(aux, stringsAsFactors = FALSE)
  aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
  colnames(aux)[1] <- "node1"
  aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)
  aux <- rbind(aux,nodes_noEdge)
  return(aux)
}


# caso de uso

pathway <- "01054"
df_pathway <- pathway2dataframe(paste0("ko", pathway))

# convertendo ko para ec pra rota ficar mais enxuta - DUVIDA

idx <- grepl("^path", df_pathway$node1) | grepl("^path", df_pathway$node2)
df_pathway <- df_pathway[!idx,]

# verificando as enzimas de um ko
# ecs <- paste0(keggLink("ec", "K01623"), collapse = ";")

df_pathway$node1 <- gsub("^ko:", "", df_pathway$node1)
df_pathway$node2 <- gsub("^ko:", "", df_pathway$node2)
#ecs <- sapply(unique(c(df_pathway$node1, df_pathway$node2)), function(i){
#  paste0(keggLink("ec", i), collapse = ";")})

#ecs[ecs==""] <- NA
#df_pathway[] <- ecs[unlist(df_pathway)]
#df_pathway <- unique(df_pathway)

uncollapse <- function(df_pathway){
  colnames(df_pathway) <- c("key", "value")
  df_pathway <- df_pathway %>%
    mutate(key = strsplit(key, ";")) %>%
    unnest(key)
  return(df_pathway)
}

df_pathway <- unique(uncollapse(df_pathway))
df_pathway <- unique(uncollapse(df_pathway[, c(2,1)]))[, c(2,1)]
colnames(df_pathway) <- c("node1", "node2")

#ecsFromPathway <- keggLink("ec", paste0("map", pathway))
#ecsFromPathway <- unname(ecsFromPathway)

#idx <- (df_pathway$node1 %in% ecsFromPathway) & (df_pathway$node2 %in% ecsFromPathway)
#df_pathway <- df_pathway[idx,]

# ecs <- keggLink("ec", "map00010")
# ecs <- unname(ecs)
# em all links temos 50 enzimas listadas
# https://www.genome.jp/dbget-bin/www_bget?pathway+map00010
# length(unique(c(df$node1, df$node2)))

pepa <- graph_from_data_frame(df_pathway, TRUE)

library(RedeR)
rdp <- RedPort()
calld(rdp)
addGraph(rdp, pepa)

plot(pepa)
