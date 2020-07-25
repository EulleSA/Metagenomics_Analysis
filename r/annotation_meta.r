
columns_csv <- c('qseqid' ,'sseqid' ,'pident' ,'length','mismatch','gapopen','qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore')
df_besthits <- read.csv("C:\\Users\\eulle\\Documents\\besthits/besthits_CoSQG.m8",header = FALSE,col.names = columns_csv,sep = "\t",stringsAsFactors = FALSE)
colnames(df_besthits) <- columns_csv

df_besthits <- subset(df_besthits, select = c("qseqid","sseqid"))

# Desmembrar os dados em 3 colunas

id <- gsub("^.*\\|(.*)\\|.*$", "\\1", df_besthits$sseqid)
id <- as.vector(gsub("GI:|gi:|UniProtKB/Swiss-Prot:|UniProtKB/TrEMBL:(.*)|p$", "\\1", id))
df_besthits["ID"] <- id

organisms <- as.vector(gsub("^.*\\|.*\\|(.*)$", "\\1", df_besthits$sseqid))
df_besthits["Organisms"] <- organisms

proteins <- as.vector(gsub("(^.*)\\|.*\\|.*$", "\\1", df_besthits$sseqid))
df_besthits["Proteins"] <- proteins


# Dividir em dois subsets

df_besthits_Uniprot <- df_besthits[grepl("[A-Z]+", df_besthits$ID),]
df_besthits_gi <- df_besthits[!grepl("[A-Z]+", df_besthits$ID),]

# Concatenar os dados
df_gi_uniprot <-  read.csv("C:\\Users\\eulle\\Documents\\gi_number-uni/ROCHA2/df_gi_uniprot.csv",header=TRUE,sep = "\t",stringsAsFactors = FALSE)

if (length(df_gi_uniprot[duplicated(df_gi_uniprot$Uniprot_ID),"Uniprot_ID"]) > 0 ){
  df_gi_uniprot <- df_gi_uniprot[!duplicated(df_gi_uniprot[ , "Uniprot_ID"]),]

}else{
  print("Não tem valores duplicados")
}

df_besthits_gi <- merge(df_besthits_gi,df_gi_uniprot,by.x = "ID",by.y = "GI",all.x = TRUE)
df_besthits_gi <- merge(df_besthits_gi,df_rev,by="Uniprot_ID",all.x = TRUE) # Uso o DF_PATHWAYS ( df_rev)
df_besthits_gi[is.na(df_besthits_gi$Pathway),"Pathway"] <- "NA"
df_besthits_gi[is.na(df_besthits_gi$KO),"KO"] <- "NA"
df_besthits_gi <- df_besthits_gi[,!(names(df_besthits_gi) %in% "Uniprot_ID")]

df_besthits_Uniprot <- merge(df_besthits_Uniprot,df_rev,by.x = "ID",by.y="Uniprot_ID",all.x = TRUE)
df_besthits_Uniprot[is.na(df_besthits_Uniprot$Pathway),"Pathway"] <- "NA"
df_besthits_Uniprot[is.na(df_besthits_Uniprot$KO),"KO"] <- "NA"


df_besthits <- rbind(df_besthits_gi,df_besthits_Uniprot)

write.table(df_besthits[df_besthits$Pathway!="NA","qseqid"],"besthits_pathways-CoSQG.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

### Gerar a imagem dos organismos

library(ggplot2)

df_besthits_organisms <- rle(sort(df_besthits$Organisms))
df_besthits_organisms <- data.frame(Organisms=df_besthits_organisms$values, n=df_besthits_organisms$lengths,stringsAsFactors = FALSE)
df_besthits_organisms_10 <- head(df_besthits_organisms[with(df_besthits_organisms,order(-n)),],n=10L)
rownames(df_besthits_organisms_10)<-NULL



df_base <- ggplot(data = df_besthits_organisms_10,aes(x=df_besthits_organisms_10$n,y=reorder(as.character(df_besthits_organisms_10$Organisms),df_besthits_organisms_10$n)))
df_base + geom_bar(stat="identity", color='skyblue',fill='steelblue') + xlab("Count of reads") + ylab("Organisms") + scale_x_continuous(limits = c(0, 35000), breaks = seq.int(0, 35000, 5000)) + ggtitle("Top 10 organisms")

## Gerar a imagem das vias.

pathways <- subset(df_besthits,df_besthits$Pathway!="NA","Pathway")
row.names(pathways) <- NULL
pathways <- strsplit(as.character(pathways$Pathway), "\\;")
pathways <- unlist(pathways,use.names=FALSE)
df_besthits_pathways <- rle(sort(pathways))

df_besthits_pathways <- data.frame(Pathways=df_besthits_pathways$values, n=df_besthits_pathways$lengths,stringsAsFactors = FALSE)
rownames(df_besthits_pathways)<-NULL

# Gerar todos

df_base_pathways <- ggplot(data = df_besthits_pathways,aes(x=df_besthits_pathways$n,y=reorder(as.character(df_besthits_pathways$Pathways),df_besthits_pathways$n)))
df_base_pathways + geom_bar(stat="identity", color='black',fill='#FF9999') + xlab("Count of reads") + ylab("Pathways") + scale_x_continuous(limits = c(0, 12000), breaks = seq.int(0, 12000, 1500))

# Gerar os 10 primeiros
df_besthits_pathways_10 <- head(df_besthits_pathways[with(df_besthits_pathways,order(-n)),],n=10L)

df_base_pathways_10 <- ggplot(data = df_besthits_pathways_10,aes(x=df_besthits_pathways_10$n,y=reorder(as.character(df_besthits_pathways_10$Pathways),df_besthits_pathways_10$n)))
df_base_pathways_10 + geom_bar(stat="identity", color='black',fill='#F49D37') + xlab("Count of reads") + ylab("Pathways") + scale_x_continuous(limits = c(0, 12000), breaks = seq.int(0, 12000, 1500)) + ggtitle("Top 10 pathways")


## Gerar a imagem da familia presente


pathways <- subset(df_besthits,df_besthits$Pathway!="NA","Pathway")
row.names(pathways) <- NULL
pathways <- strsplit(as.character(pathways$Pathway), "\\;")
pathways <- unlist(pathways,use.names=FALSE)


df_biosurfdb_pathways <- read.csv("C:\\Users\\eulle\\OneDrive/Área de Trabalho/Biosurfdb.csv",header=TRUE,sep = ",",stringsAsFactors = FALSE)
pathways <- data.frame(Pathways=pathways,stringsAsFactors = FALSE) 
pathways <- merge(pathways,df_biosurfdb_pathways,by.x="Pathways",by.y="name",all.x = TRUE)
pathways <- pathways[!is.na(pathways$family),"family"]
pathways <- rle(sort(pathways))
family_counts <- data.frame(Family=pathways$values,counts=pathways$lengths)

family_counts_surfc_degrad <- subset(family_counts,Family=="Degradation" | Family=="Surfactants")

df_base_family <- ggplot(data = family_counts_surfc_degrad,aes(x=family_counts_surfc_degrad$counts,y=reorder(as.character(family_counts_surfc_degrad$Family),family_counts_surfc_degrad$counts)))
df_base_family + geom_bar(stat="identity", color='black',fill='#FF9999') +ggtitle("Pathways grouped by Biosurfdb's families") + xlab("Count of Pathways") + ylab("Biosurfdb Families") + scale_x_continuous(limits = c(0, 25000), breaks = seq.int(0, 25000, 5000))

# GERAR GRÁFICO SOBRE O QUÃO COMPLETADO ESTÁ A ROTA

df_kos_count_plot <- ggplot(data=df_kos_count,aes(x=df_kos_count$KO_per_total,y=reorder(df_kos_count$Name,df_kos_count$KO_per_total)))
df_kos_count_plot + geom_bar(stat="identity", color='black',fill='#FF9999') + xlab("Count of KO in percentage") + ylab("Pathways")+ggtitle("Count of pathways per compunds") + scale_x_continuous(limits = c(0, 30), breaks = seq.int(0, 30, 10))


## GERAR A IMAGEM DAS CONTAGENS 10 primeiros
df_kos_count_10 <- head(df_kos_count[with(df_kos_count,order(-KO_per_total)),],n=10L)
df_kos_count_plot <- ggplot(data=df_kos_count_10,aes(x=df_kos_count_10$KO_per_total,y=reorder(df_kos_count_10$Name,df_kos_count_10$KO_per_total)))
df_kos_count_plot + geom_bar(stat="identity", color='black',fill='#FF9999') + xlab("KO rate") + ylab("Pathways") + scale_x_continuous(limits = c(0, 30), breaks = seq.int(0, 30, 10)) + ggtitle("KO rate (found/pathway total)")




