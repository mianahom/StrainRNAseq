library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")
library("ggplot2")
library("clusterProfiler")
library("enrichplot")
library("ggupset")

#count and metadata

setwd("/Users/mianahom/Downloads/Strain_RNAseq/FemaleResults")
directory <- "/Users/mianahom/Downloads/Strain_RNAseq/counts"
list.files(directory)
sampleFiles <- list.files(directory, pattern="*.counts$")

meta <- read.csv("../meta.csv")
meta
all(sampleFiles == meta$Files)

sampleTable <- data.frame(
  sampleName = meta$Sample,
  fileName = sampleFiles,
  sex = meta$Sex,
  condition = meta$condition
)

#MaleTable <- subset(sampleTable, sex=="Male")

FemaleTable <- subset(sampleTable, sex=="Female")


ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = FemaleTable, 
  directory = directory, 
  design = ~ condition
)

ddsHTSeq$condition

# To replace the order with one of your choosing, create a vector with the order you want:
conditions <- c("wild","heterozygous")

# Then reset the factor levels:
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = conditions)

# verify the order
ddsHTSeq$condition

# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 20
keep <- sumcounts > 50

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)
resultsNames(dds)
res <- results(dds, name="condition_heterozygous_vs_wild")
summary(res)

res_shrink <- lfcShrink(dds,type="ashr",coef="condition_heterozygous_vs_wild")

res_ranked <- as.data.frame(res_shrink) %>%
  filter(padj < 0.1) %>%
  arrange(-abs(log2FoldChange))

#mart /dataset
# ensembl host:
# to list archived version hosts: listEnsemblArchives()
listEnsemblArchives()
ensemblhost <- "https://oct2024.archive.ensembl.org"

listMarts(host=ensemblhost)

# create an object for the Ensembl Genes v113 mart
mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", mirror = 'asia' )

# occasionally ensembl will have connectivity issues. we can try an alternative function:
# select a mirror: 'www', 'uswest', 'useast', 'asia'
# mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

# see a list of datasets within the mart

listDatasets(mart)

# figure out which dataset is mouse
# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="mmusculus_gene_ensembl")

# there's only one match, get the name
mousedata <- searchDatasets(mart,pattern="mmusculus_gene_ensembl")[,1]

# create an object for the mouse dataset
mouse_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror='asia', dataset = mousedata)

# if above there were connectivity issues and you used the alternative function then:
# select a mirror: 'www', 'uswest', 'useast', 'asia'
# killi_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = killidata, mirror = "useast")

#########################
# Query the mart/dataset
#########################

# filters, attributes and values

# see a list of all "filters" available for the mummichog dataset.
# at the time of writing, over 300
listFilters(mouse_mart)

# see a list of all "attributes" available

listAttributes(mart = mouse_mart, page="feature_page")

# we can also search the attributes and filters
searchAttributes(mart = mouse_mart, pattern = "ensembl_gene_id")

searchFilters(mart = mouse_mart, pattern="ensembl")

# get gene names and transcript lengths when they exist
annotation <- 
  getBM(
    filter="ensembl_gene_id",
    value=rownames(res),
    attributes=c("ensembl_gene_id","description","transcript_length","external_gene_name"),
    mart=mouse_mart
  )

# pick only the longest transcript for each gene ID
annotation <- 
  group_by(annotation, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),external_gene_name=unique(external_gene_name),transcript_length=max(transcript_length)) %>%
  as.data.frame()

# get GO term info
# each row is a single gene ID to GO ID mapping, so the table has many more rows than genes in the analysis
go_annotation <- getBM(
  filter="ensembl_gene_id",
  value=rownames(res),
  attributes=c("ensembl_gene_id","description","go_id","name_1006","definition_1006","namespace_1003"),
  mart=mouse_mart
)

#let's look at the annotation:
filter(annotation, description!="") %>% head()
head(annotation)
# and Go term:
head(go_annotation)
# put results and annotation in the same table
# (Add gene names to DE results)

res_shrink_data <- res_shrink %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_gene_id") 

res_shrink_ann <- merge(x=res_shrink_data,y=annotation, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x=TRUE)


res_data <- res %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_gene_id") 

res_ann <- merge(x=res_data,y=annotation, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x=TRUE)



##Plots to represent the data
#dispersion estimates
png("DispersionEstimates.png", width = 1200, height = 1200, res = 300)
plotDispEsts(dds)
dev.off()

#l2fc vs shrunken l2fc
png("L2FC_vsShrunkL2FC.png", width = 1200, height = 1200, res = 300)
data.frame(l2fc=res$log2FoldChange, l2fc_shrink=res_shrink$log2FoldChange, padj=res$padj) %>%
  filter(l2fc > -5 & l2fc < 5 & l2fc_shrink > -5 & l2fc_shrink < 5) %>%
  ggplot(aes(x=l2fc, y=l2fc_shrink,color=padj < 0.1)) +
  geom_point(size=.25) + 
  geom_abline(intercept=0,slope=1, color="gray")
dev.off()

# MA plot
png("MAplot.png", width = 1200, height = 1200, res = 300)
plotMA(res, ylim=c(-4,4))
dev.off()

png("MAplot_shrunk.png", width = 1200, height = 1200, res = 300)
plotMA(res_shrink, ylim=c(-4,4))
dev.off()
##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.5,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)


png("PCA.png", width = 1200, height = 1200, res = 300)
plotPCA(vsd, intgroup="condition")
dev.off()

#QC check target gene
plotCounts(dds,gene = "ENSMUSG00000031393")



##HEATMAPS
dataframe <- data.frame(colData(dds)[,c("condition")])
rownames(dataframe) <- colnames(dds)
colnames(dataframe) <- c("condition")

# pull the top 50 genes by shrunken log2 fold change
rs <- res_shrink_ann 

rs <- rs %>%
  mutate(external_gene_name = ifelse(external_gene_name == "", ensembl_gene_id, external_gene_name))


top50 <- rs %>% 
  data.frame() %>%
  filter(padj < 0.1) %>%
  arrange(-abs(log2FoldChange)) %>%
  head(n = 50) %>%
  select("ensembl_gene_id","external_gene_name")


# heatmap of normalized variance-stabilized transformed counts top 50
png("top50_heatmap.png", width = 5000, height = 5000, res = 500)
pheatmap(
  assay(vsd)[top50$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=top50$external_gene_name,
  annotation_col=dataframe
)
dev.off()

#spearman corr top 50
png("top50_spearmancorrelation_heatmap.png", width = 5000, height = 5000, res = 500)
pheatmap(
  assay(vsd)[top50$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=top50$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()
#rescaled by overall expression level
rescaled <- assay(vsd)-rowMeans(assay(vsd))

#rescaled top 50 spearman corr
png("top50rescaled_spearmancorrelation_heatmap.png", width = 5000, height = 5000, res = 500)
pheatmap(
  rescaled[top50$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=top50$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

candidates <- c("Kcnk3","Kcnk9","Th",
                "Gfap","Phox2b","Epas1",
                "Hif1a","Fos","Egln1",
                "Kcna3","Kcnd3","Hcn1",
                "Cacna1c","Scn10a","Scn11a",
                "Ndufs2","Cox4i2","Cox8b",
                "Hmox1","Ddc","Drd2",
                "P2rx2","P2rx3",
                "Lepr", "Trpm7","Or51e2")


candidates <- data.frame(values=candidates)
names(candidates)[1] <- "external_gene_name"
candidates <- left_join(candidates, annotation, by = "external_gene_name")

Hypoxia <- c("Epas1", "Hif1a", "Egln1", "Fos") %>%
  data.frame(external_gene_name = .) %>%  # Create a data frame with column name
  left_join(annotation, by = "external_gene_name")

IonElectro <- c("Kcnk3","Kcnk9","Kcna3","Kcnd3","Hcn1",
                "Cacna1c","Scn10a","Scn11a") %>%
  data.frame(external_gene_name = .) %>%  
  left_join(annotation, by = "external_gene_name")

Mitochondrial <- c("Ndufs2","Cox4i2","Cox8b","Hmox1") %>%
  data.frame(external_gene_name = .) %>%  
  left_join(annotation, by = "external_gene_name")

Neurotransmitter <- c("Th","Ddc","Drd2","P2rx2","P2rx3",
                      "Lepr", "Trpm7","Or51e2") %>%
  data.frame(external_gene_name = .) %>%  
  left_join(annotation, by = "external_gene_name")

# heatmap of normalized variance=stabilized tranformed counts 


png("candidate_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  assay(vsd)[candidates$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=candidates$external_gene_name,
  annotation_col=dataframe
)
dev.off()

#spearman corr top 50
png("candidate_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  assay(vsd)[candidates$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=candidates$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

#rescaled top 50 spearman corr
png("candidaterescaled_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  rescaled[candidates$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=candidates$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

png("Hypoxiarescaled_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  rescaled[Hypoxia$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=Hypoxia$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

png("IonElectrorescaled_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  rescaled[IonElectro$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=IonElectro$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

png("Mitochondrialrescaled_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  rescaled[Mitochondrial$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=Mitochondrial$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()

png("Neurotransmitterrescaled_spearman_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  rescaled[Neurotransmitter$ensembl_gene_id,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  labels_row=Neurotransmitter$external_gene_name,
  annotation_col=dataframe,
  clustering_distance_rows="correlation"
)
dev.off()


#enrichment analyses:
# get ENSEMBL gene IDs for genes with padj < 0.1
genes <- rownames(res[which(res$padj < 0.1),])
# get ENSEMBL gene IDs for universe (all genes with non-NA padj passed independent filtering)
univ <- rownames(res[!is.na(res$padj),])
# pull out the columns of go_ann containing GO IDs and descriptions, keep only unique entries. 
gonames <- unique(go_annotation[,c(3,4)]) %>% unique()

enrich <- enricher(
  gene=genes,
  universe=univ,
  TERM2GENE=go_annotation[,c(3,1)],
  TERM2NAME=gonames
)

enrichobj<-as.data.frame(enrich)
# extract log2 fold changes, only for genes that passed independent filtering. 
l2fcs <- as.data.frame(res) %>% 
  filter(!is.na(padj)) %>%
  select(log2FoldChange)

# put log2FCs in a vector, add gene IDs as names, sort 
l2fcvec <- l2fcs[,1]
names(l2fcvec) <- rownames(l2fcs)
l2fcvec <- sort(l2fcvec, decreasing=TRUE)

res_gsea <- GSEA(
  geneList=l2fcvec,
  TERM2GENE=go_annotation[,c(3,1)],
  TERM2NAME=gonames
) 

gseaplot(res_gsea, by = "all", title = enrich$Description[2], geneSetID = 2)

gseaplot(res_gsea, by = "all", title = enrich$Description[1], geneSetID = 1)

gseaobj <-as.data.frame(res_gsea)


png("OverRepresentation_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(enrich)
dev.off()

png("OverRepresentation_barplot.png", width = 3000, height = 2000, res = 300)
barplot(enrich)
dev.off()

png("GSEA_dotplot.png", width = 2000, height = 2000, res = 300)
dotplot(res_gsea, showCategory = 10) +
  ggtitle("GSEA Dot Plot")
dev.off()


# Example data for heatmap
#pathwayGenes <- pathway[["Example_Set"]]
#rankedGenes <- ranks[pathwayGenes]

# Create a heatmap
#pheatmap(matrix(rankedGenes, nrow = 1), cluster_cols = FALSE,
#    labels_col = names(rankedGenes),
#    main = "Gene Set Heatmap")
png("OverRepresentation_upsetplot.png", width = 5000, height = 2000, res = 300)
upsetplot(enrich)
dev.off()

png("GSEA_upsetplot.png", width = 5000, height = 2000, res = 300)
upsetplot(res_gsea)
dev.off()

# truncate GO descriptions if you want (for cleaner plotting lataer)
#gonames[,2] <- sapply(substring(gonames[,2],1,30), FUN=function(x){if(nchar(x) == 30){x <- paste(x,"...",sep="")}; x})

png("GSEA_ridgeplot.png", width = 10000, height = 5000, res = 500)
ridgeplot(res_gsea,label_format=100)
dev.off()


#counts for regular gene names
df <- as.data.frame(counts(dds,normalized =TRUE))
mousegenes <- annotation %>% select(ensembl_gene_id, external_gene_name)
df <- merge(x=mousegenes,y=df, by.x="ensembl_gene_id", by.y=0, all.y=TRUE)
colnames(df)[1] <- "ens_gene_id"
ft <- data.frame(treatment=dds$condition)


countPlots <- function(df, ens_gene_id, factor_table){
  # row in df containing gene of interest
  gene <- which(df[,1]==ens_gene_id)
  # get counts
  counts <- df[gene, 3:dim(df)[2]] %>% unlist() %>% as.numeric()
  counts <- log(counts+1, base=10)
  
  # df containing counts and factors
  newdf <- data.frame(count_data=counts, factor_table)
  
  ggplot(newdf, aes(x=newdf[,2],y=newdf[,1])) +
    geom_boxplot() + 
    geom_jitter() +
    ggtitle(df[gene,]$external_gene_name) +
    ylab("log-scaled (log 10 +1 ) normalized counts") +
    xlab(colnames(newdf)[2])
}

#MeCP2
#ENSMUSG00000031393
png("male_countplot_Mecp2.png", width = 1000, height = 1000, res = 300)
countPlots(df, "ENSMUSG00000031393", ft)
dev.off()

for (i in 1:26) { 
  name <- candidates$external_gene_name[i] 
  file_name <- paste0("male_countplot_",name,".png")
  png(file_name, width = 1000, height = 1000, res = 300)
  countplot <- countPlots(df, candidates$ensembl_gene_id[i], ft)
  print(countplot)
  dev.off()
}

#Save results:
#normalized counts and gene names:
write.csv(df,file="MaleNormalizedCounts.csv")
#results and annotation
write.csv(res_ann, file="MaleResults.csv")
write.csv(go_annotation, file="GoTerms.csv")
#GSEA and Overrepresentation
write.csv(enrichobj, file="OverRepresentation.csv")
write.csv(gseaobj, file="GSEA.csv")
write.csv(univ, file="backgroundgenes.csv")
write.csv(genes, file="GenesforShiny.csv")

clipr::write_clip(univ)
clipr::write_clip(genes)