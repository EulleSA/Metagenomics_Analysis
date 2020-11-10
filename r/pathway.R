#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("pathview")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("XML",type="binary")
library(pathview)
library(AnnotationDbi)
library(stringr)


for(index in 1:length(df_kos_count$cod_pathway)){
  
  pathway <- as.character(df_kos_count[[2]][index]) 
  IDs <- str_split(df_kos_count[df_kos_count$cod_pathway==pathway,][4],";")
  IDs <- unlist(IDs, use.names=FALSE)
  #now <- format(Sys.time(), "%Y-%m-%d--%H-%M-%S")
  img <- pathview(gene.data = IDs,
                  pathway.id = pathway,
                  out.suffix = df_kos_count[[1]][index],
                  species = "ko",
                  high = list(gene = "darkseagreen1"),
                  kegg.native = TRUE,
                  same.layer = FALSE,
                  new.signature = FALSE,
                  plot.col.key = FALSE,
                  map.symbol = TRUE,
                  gene.annotpkg = NA,
                  map.null = FALSE)
}

junk<- dir(path="./",pattern = "*.xml|.$")
file.remove(junk)
