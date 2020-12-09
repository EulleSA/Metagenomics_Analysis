devtools::install_github("daniloimparato/easylayout", ref = "dadamorais")

library(igraph)
library(ggraph)
library(scatterpie)
library(easylayout)
library(stringr)
#load("grafo.Rdata")


# Criar um dataframe
complete_dataframe <- function(rowname_taxonomy,cod_pathway){
  
  #df_kos_count <- 
  #rowname_taxonomy <- df_taxonomybesthits[unique(unlist(str_split(df_taxonomybesthits$KEGG_ko,","))) %in% str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]] ,c("KEGG_ko","best_tax_level")]
  rowname_taxonomy <- df_taxonomybesthits[unlist(lapply(str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]], function(x) grep(x, df_taxonomybesthits$KEGG_ko, fixed = TRUE))),c("KEGG_ko","best_tax_level")]
  
  
  name_org <- unique(rowname_taxonomy$best_tax_level)
  colname_kos <- V(pepa)$name
  kos_binary <- lapply(name_org,function(name_org){as.numeric(colname_kos %in% rowname_taxonomy[rowname_taxonomy$best_tax_level==name_org,"KEGG_ko"])})
  
  if(length(kos_binary) == 0){
    df_categories <- data.frame(matrix(ncol = length(colname_kos), nrow = 0))
    colnames(df_categories) <- colname_kos
    
    `%notin%`<- Negate(`%in%`)
    unclassified_kos <- str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]]
    unclassified_kos <- Unclassified_kos[unclassified_kos %notin% rowname_taxonomy$KEGG_ko]
    
    df_categories <- rbind(df_categories,ifelse(colname_kos %in% unclassified_kos,1,0))
    rownames(df_categories)[nrow(df_categories)] <- "Unclassified"
    
    df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>=1,0,1))
    rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  }
  else{
    df_categories <- data.frame(matrix(unlist(kos_binary), nrow=length(kos_binary), byrow=T),stringsAsFactors=FALSE)
    
    `%notin%`<- Negate(`%in%`)
    unclassified_kos <- str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]]
    unclassified_kos <- unclassified_kos[unclassified_kos %notin% rowname_taxonomy$KEGG_ko]
    
    rownames(df_categories) <- name_org
    colnames(df_categories) <- V(pepa)$name
    
    df_categories <- rbind(df_categories,ifelse(colnames(df_categories) %in% unclassified_kos,1,0))
    rownames(df_categories)[nrow(df_categories)] <- "Unclassified"
    
    df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>=1,0,1))
    rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  }
  # df_categories <- data.frame(matrix(unlist(kos_binary), nrow=length(kos_binary), byrow=T),stringsAsFactors=FALSE)
  # 
  # `%notin%`<- Negate(`%in%`)
  # Unclassified_kos <- str_split(df_kos_count[df_kos_count$cod_pathway == cod_pathway,"KO_found"],";")[[1]]
  # Unclassified_kos <- Unclassified_kos[Unclassified_kos %notin% rowname_taxonomy$KO]
  # 
  # rownames(df_categories) <- name_org
  # colnames(df_categories) <- V(pepa)$name
  # 
  # df_categories <- rbind(df_categories,ifelse(colnames(df_categories) %in% Unclassified_kos,1,0))
  # rownames(df_categories)[nrow(df_categories)] <- "Unclassified"
  # 
  # df_categories <- rbind(df_categories,ifelse(colSums(df_categories)>=1,0,1))
  # rownames(df_categories)[nrow(df_categories)] <- "Not Found"
  
  
  
  return(df_categories)
}

df_taxonomybesthits <- read.csv("D:/Downloads/taxmy_CoSqG",sep = "\t",stringsAsFactors = FALSE)


# names(df_taxonomybesthits)[names(df_taxonomybesthits) %in% c("best_tax_level","X.query_name")] <- c("query_name","taxonomy")
# 
# df_taxonomybesthits <- subset(df_taxonomybesthits,KEGG_Pathway != "",select = c("query_name","taxonomy","KEGG_ko","KEGG_Pathway"))
# df_taxonomybesthits$KEGG_ko <- gsub("ko:","",df_taxonomybesthits$KEGG_ko)
# 
# write.table(df_taxonomybesthits,"D:/Downloads/taxbesthitsby26.csv",row.names = FALSE,col.names = TRUE)

rowname_taxonomy <- df_taxonomybesthits[unlist(lapply(str_split(df_kos_count[df_kos_count$cod_pathway == "ko00281","KO_found"],";")[[1]], function(x) grep(x, df_taxonomybesthits$KEGG_ko, fixed = TRUE))),c("KEGG_ko","best_tax_level")]

# Gerando layout

layout <- easylayout::vivagraph(pepa)


df_organisms <- complete_dataframe(df_taxonomybesthits,"ko00281")

# Inserindo dados nos vertices
# Categoria é um valor binário onde 1 - COR e 0 - NÃO COR
# x e y são as coordenadas
# Ver esses negócios dA LISTA
categories <- list(name=V(pepa)$name)

categories <- append(categories,lapply(lapply(setNames(split(df_organisms,seq(nrow(df_organisms))),rownames(df_organisms)),unname),unlist))
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
    ,arrow = arrow(length = unit(1.5, 'mm'), type = "closed")
    ,alpha = 0.5
  ) +
  
  # Pie charts
  geom_scatterpie(
    data    = igraph::as_data_frame(pepa, "vertices")
    ,mapping = aes(x, y, group = name, r = 7)
    ,cols    = rownames(df_organisms)
    ,color   = "black"
  ) +
  
  # Aspect ratio
  coord_fixed() +
  # Node label
  geom_node_text(mapping = aes(x, y), label = gsub("^ec:", "", V(pepa)$name),
                 nudge_y = 15) +
  theme_void()


