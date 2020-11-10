devtools::install_github("daniloimparato/easylayout", ref = "dadamorais")

library(igraph)
library(ggraph)
library(scatterpie)
library(easylayout)
library(stringr)
#load("grafo.Rdata")

# Gerando layout
layout <- easylayout::vivagraph(pepa)

# Criar um dataframe
complete_dataframe <- function(rowname_taxonomy,cod_pathway){
  
  rowname_taxonomy <- df_taxonomybesthits[df_taxonomybesthits$KO %in% str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]] ,c("KO","taxonomy")]
  name_org <- unique(rowname_taxonomy$taxonomy)
  colname_kos <- V(pepa)$name
  kos_binary <- lapply(name_org,function(name_org){as.numeric(colname_kos %in% rowname_taxonomy[rowname_taxonomy$taxonomy==name_org,"KO"])})
  df_categories <- data.frame(matrix(unlist(kos_binary), nrow=length(kos_binary), byrow=T),stringsAsFactors=FALSE)
  rownames(df_categories) <- name_org
  colnames(df_categories) <- V(pepa)$name
  df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>1,0,1))
  rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  
  return(df_categories)
}


df_taxonomybesthits <- read.csv("C:/Users/eulle/Documents/metagenomica/taxonomy_besthits-CoSQG.csv",sep = "]")
df_categories <- complete_dataframe(df_taxonomybesthits,"ko01054")


# rowname_taxonomy <- df_taxonomybesthits[df_taxonomybesthits$KO %in% str_split(df_kos_count[df_kos_count$cod_pathway == "ko00071","KO_found"],";")[[1]] ,c("KO","taxonomy")]
# name_org <- unique(rowname_taxonomy$taxonomy)
# colname_kos <- V(pepa)$name
# kos_binary <- lapply(name_org,function(name_org){as.numeric(colname_kos %in% rowname_taxonomy[rowname_taxonomy$taxonomy==name_org,"KO"])})
# df_categories <- data.frame(matrix(unlist(kos_binary), nrow=length(kos_binary), byrow=T),stringsAsFactors=FALSE)
# rownames(df_categories) <- name_org
# colnames(df_categories) <- V(pepa)$name
# df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>1,0,1))

#rowname_taxonomy <- df_taxonomybesthits[df_taxonomybesthits$KO %in% str_split(df_kos_count[df_kos_count$cod_pathway == "ko00071","KO_found"],";")[[1]] ,c("KO","taxonomy")]
#name_org <- unique(rowname_taxonomy$taxonomy)


#df_orgbyko <- data.frame(row.names = unique(rowname_taxonomy$taxonomy),stringsAsFactors = FALSE)
#colnames(df_orgbyko)<- colname_kos
         


# Inserindo dados nos vertices
# Categoria é um valor binário onde 1 - COR e 0 - NÃO COR
# x e y são as coordenadas
# Ver esses negócios dA LISTA
categories <- list(name=V(pepa)$name)

categories <- append(categories,lapply(lapply(setNames(split(df_categories,seq(nrow(df_categories))),rownames(df_categories)),unname),unlist))
categories[["x"]] <- layout[,1]
categories[["y"]] <- layout[,2]
vertex_attr(pepa) <- categories

# Inserindo o tamanho dos nós
# PERGUNTAR COMO COLOCA O MESMO TAMANHO PARA TODOS 
#V(pepa)$size <- V(pepa)$categoria 

# Plotando
ggraph(pepa, x = x, y = y) +
  
  # Arestas primeiro
  geom_edge_link(
    aes(end_cap = circle(1.5 + 0.5, 'mm'))
    ,color = "#808080"
    ,arrow = arrow(length = unit(1, 'mm'), type = "closed")
    ,alpha = 0.5
  ) +
  
  # Pie charts
  geom_scatterpie(
    data    = igraph::as_data_frame(pepa, "vertices")
    ,mapping = aes(x, y, group = name, r = 5)
    ,cols    = rownames(df_categories)
    ,color   = "black"
  ) +
  
  # Aspect ratio
  coord_fixed() +
  # Node label
  geom_node_text(mapping = aes(x, y), label = gsub("^ec:", "", V(pepa)$name),
                 nudge_y = 15) +
  theme_void()


