
## Gene-level differential expression analysis using DESeq2

rm(list=ls())
setwd("~/Desktop/Adrenal/")

## Load required libraries


library("ggplot2")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggbeeswarm")
library("reshape2")
library("pheatmap")
library("genefilter")
library("clusterProfiler")
library("org.Hs.eg.db")
library("Homo.sapiens")
library("limma")
library("sva")
library("ggthemes")
library("ggdendro")
library("rmarkdown")
library("ggpubr")
library("tidyverse")
library("DEGreport")
library("factoextra")
library('fgsea')
library('DT')
library('rgl')
library("biomaRt")

## Load data


phenodata <- read.table("metadata/phenodata.csv", sep=",", header=TRUE)

#phenodata <- phenodata %>% filter (!(Number == 9)) Removing Liver sample DOES NOT change PC2-3 dramatically
View(phenodata)


countdata <- read.table(file="data/counts.csv", header=TRUE, sep=",", row.names=1)

#countdata <- countdata %>% select (-L13514) Remove liver sample

View(countdata)

ncol(countdata) == nrow(phenodata)

colnames(countdata) ==  phenodata$Sample # They should be in the same order


# Convert to matrix
class(countdata)
countdata <- as.matrix(countdata)

# Assign condition according to colnames
tissue <- phenodata$Tissue
tissue

# It is better in R if the first level of a factor is the reference level:
tissue %<>% relevel("Control")
tissue

batch <- factor(phenodata$Batch)
batch

stage <- factor(phenodata$Stage)
stage

sex <- factor(phenodata$Sex)
sex %<>% relevel("XY")
sex

origin <- factor(phenodata$Origin)
origin

sample.stage <- phenodata$Sample_stage
sample.stage

phenodata$sample.dev.origin <- paste(phenodata$Stage, phenodata$Sample, sep = "_")
sample.dev.origin <- phenodata$sample.dev.origin
rownames(phenodata) <- phenodata$sample.dev.origin
View(phenodata)

number <- phenodata$Number
number


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), batch, tissue, sex, stage, origin, number, sample.dev.origin)
coldata



dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~batch + sex + stage + tissue)

dds

dds <- dds[ rowSums(counts(dds)) > 1, ] #We apply the most minimal filtering rule: removing rows of the DESeqDataSet that have no counts, or only a single count across all samples. 
nrow(dds)

dds <- DESeq(dds)

sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

head(normalized_counts)

write.table(normalized_counts, file="results/Adrenal-Controls/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# COUNT DATA QC : VST Transformation and PCA Plots #

# vst is faster to compute and is less sensitive to high count outliers. Is recommended for large datasets.

vsd <- vst(dds, blind=TRUE)
head(assay(vsd),1)
hist(assay(vsd))
names(colData(dds))


plotPCA.tissue <- plotPCA(vsd, intgroup=c("tissue"),ntop = Inf) + theme_classic() 
plotPCA.stage <- plotPCA(vsd, intgroup=c("stage"),ntop = Inf)+ theme_classic()
plotPCA.batch <- plotPCA(vsd, intgroup=c("batch"),ntop = Inf)+ theme_classic()
plotPCA.sex <- plotPCA(vsd, intgroup=c("sex"),ntop = Inf)+ theme_classic()
plotPCA.origin <- plotPCA(vsd, intgroup=c("origin"),ntop = Inf) + theme_classic()

pcaData.tissue <- plotPCA(vsd, intgroup=c("tissue", "stage"), returnData=TRUE, ntop = Inf)
percentVar <- round(100 * attr(pcaData.tissue, "percentVar"))
plot.tissue.stage <- ggplot(pcaData.tissue, aes(PC1, PC2, color=tissue, shape=stage)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text(aes(label=origin), hjust=-0.25, vjust=-0.25, size=3) + theme_classic()

### Plot for rest of components #######

vsd_mat <- assay(vsd)
pca.all <- prcomp(t(vsd_mat))
Screeplot.all <- fviz_eig(pca.all)

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(phenodata, pca.all$x)
df

ggplot(df) + geom_point(aes(x=PC1, y=PC3, color = Stage), size = 3)+ theme_classic()
ggplot(df) + geom_point(aes(x=PC2, y=PC3, color = Stage), size = 3)+ theme_classic()

ggplot(df) + geom_point(aes(x=PC1, y=PC3, color = Sex), size = 3)+ theme_classic()
ggplot(df) + geom_point(aes(x=PC2, y=PC3, color = Sex), size = 3)+ theme_classic()

ggplot(df) + geom_point(aes(x=PC1, y=PC3, color = Origin), size = 3)+ theme_classic()
plotPCA.PC2.PC3 <- ggplot(df) + geom_point(aes(x=PC2, y=PC3, color = Origin), size = 3)+ theme_classic()

ggplot(df, aes(x=PC1, y=PC3, color = batch), size = 3) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
  coord_fixed()+ geom_text(aes(label=origin), hjust=-0.25, vjust=-0.25, size=3) + theme_classic()

ggplot(df, aes(x=PC2, y=PC3, color = batch), size = 3) + geom_point(size=3) + xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) + 
  coord_fixed()+ geom_text(aes(label=origin), hjust=-0.25, vjust=-0.25, size=3) + theme_classic()


###### only Adrenal samples ###########

stages.adrenal <- phenodata %>% filter (Tissue == 'Adrenal')

View(stages.adrenal)

col_order.adrenal <- as.vector(stages.adrenal$Sample)

countdata.adrenal <- countdata[, col_order.adrenal]

col_order.adrenal %in% colnames(countdata.adrenal)

col_order.adrenal == colnames(countdata.adrenal)

colnames(countdata.adrenal)

# Assign condition according to colnames
names(stages.adrenal)

stages.adrenal$stage.sex <- factor(paste0(stages.adrenal$Stage, stages.adrenal$Sex))
stage.sex <- stages.adrenal$stage.sex

batch.adrenal <- factor(stages.adrenal$Batch)
batch.adrenal

stage.adrenal <- factor(stages.adrenal$Stage)
stage.adrenal

sex.adrenal <- factor(stages.adrenal$Sex)
sex.adrenal %<>% relevel("XY")
sex.adrenal

sample.adrenal.dev.origin <- stages.adrenal$sample.dev.origin

rownames(stages.adrenal) <- stages.adrenal$Sample
stages.adrenal

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata.adrenal <- data.frame(row.names=colnames(countdata.adrenal), batch.adrenal, sex.adrenal, stage.adrenal, sample.adrenal.dev.origin, stage.sex)

coldata.adrenal



dds.adrenal <- DESeqDataSetFromMatrix(countData=countdata.adrenal, colData=coldata.adrenal, design=~batch.adrenal + sex.adrenal + stage.adrenal)

dds.adrenal

dds.adrenal <- dds.adrenal[ rowSums(counts(dds.adrenal)) > 1, ] #We apply the most minimal filtering rule: removing rows of the DESeqDataSet that have no counts, or only a single count across all samples. 
nrow(dds.adrenal)


dds.adrenal <- DESeq(dds.adrenal)
sizeFactors(dds.adrenal)
normalized_counts.adrenal <- counts(dds.adrenal, normalized=TRUE)
head(normalized_counts.adrenal)

write.table(normalized_counts.adrenal, file="results/Adrenal/normalized_counts.adrenal.txt", sep="\t", quote=F, col.names=NA)


# COUNT DATA QC : VST Transformation and PCA Plots #

vsd.adrenal <- vst(dds.adrenal, blind=TRUE)
head(assay(vsd.adrenal),1)
hist(assay(vsd.adrenal))
names(colData(dds.adrenal))

plotPCA.stage.adrenal <- plotPCA(vsd.adrenal, intgroup=c("stage.adrenal"), ntop = Inf)+ theme_classic()
plotPCA.batch.adrenal <- plotPCA(vsd.adrenal, intgroup=c("batch.adrenal"), ntop = Inf)+ theme_classic()
plotPCA.sex.adrenal <- plotPCA(vsd.adrenal, intgroup=c("sex.adrenal"), ntop = Inf)+ theme_classic()


pcaData.tissue.adrenal <- plotPCA(vsd.adrenal, intgroup=c("batch.adrenal", "sex.adrenal", "sample.adrenal.dev.origin"), ntop = Inf, returnData=TRUE)
percentVar.adrenal <- round(100 * attr(pcaData.tissue.adrenal, "percentVar"))
plot.tissue.stage.adrenal <- ggplot(pcaData.tissue.adrenal, aes(PC1, PC2, color=stage.adrenal, shape=batch.adrenal)) +
  geom_point(size=3) + theme_classic() +
  xlab(paste0("PC1: ",percentVar.adrenal[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.adrenal[2],"% variance")) + 
  coord_fixed()

plot.stage.adrenal <- ggplot(pcaData.tissue.adrenal, aes(PC1, PC2, color=stage.adrenal)) +
  geom_point(size=3) + theme_classic() +
  xlab(paste0("PC1: ",percentVar.adrenal[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.adrenal[2],"% variance")) 

plot.sex.stage.adrenal <- ggplot(pcaData.tissue.adrenal, aes(PC1, PC2, color=stage.adrenal, shape=sex.adrenal)) +
  geom_point(size=3) + theme_classic() +
  xlab(paste0("PC1: ",percentVar.adrenal[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.adrenal[2],"% variance")) 

plot.sex.stage.adrenal.labels <- ggplot(pcaData.tissue.adrenal, aes(PC1, PC2, color=stage.adrenal, shape=sex.adrenal)) +
  geom_point(size=3) + theme_classic() +
  xlab(paste0("PC1: ",percentVar.adrenal[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.adrenal[2],"% variance")) + geom_text(aes(label=stages.adrenal$Sample), hjust=-0.25, vjust=0, size=2) 

### Plot for rest of components #######

vsd_mat.adrenal <- assay(vsd.adrenal)
pca.adrenal <- prcomp(t(vsd_mat.adrenal))
Screeplot.adrenal <- fviz_eig(pca.adrenal)

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df.adrenal <- cbind(stages.adrenal, pca.adrenal$x)

pca.adrenal.2.3.labels <- ggplot(df.adrenal, aes(PC2, PC3, color = Sex, shape = stage.adrenal)) + geom_point(size=3)+ 
  coord_fixed() + theme_classic() + geom_text(aes(label=stages.adrenal$Sample), hjust=-0.25, vjust=0, size=3) 

pca.adrenal.2.3 <- ggplot(df.adrenal, aes(PC2, PC3, color = Sex, shape = stage.adrenal)) + geom_point(size=3)+ 
  coord_fixed() + theme_classic()

pca.adrenal.2.3.labels.batch <- ggplot(df.adrenal, aes(PC2, PC3, color = as.factor(stages.adrenal$Batch), shape = stage.adrenal)) + geom_point(size=3)+ coord_fixed() + theme_classic() + geom_text(aes(label=stages.adrenal$Sample), hjust=-0.25, vjust=0, size=3) 

pca.adrenal.2.3.batch <- ggplot(df.adrenal, aes(PC2, PC3, color = as.factor(stages.adrenal$Batch), shape = stage.adrenal)) + geom_point(size=3)+ 
  coord_fixed() + theme_classic()

ggplot(df.adrenal) + geom_point(aes(x=PC1, y=PC3, color = as.factor(stages.adrenal$Batch)), size = 3)+ theme_classic()
ggplot(df.adrenal) + geom_point(aes(x=PC2, y=PC3, color = as.factor(stages.adrenal$Batch)), size = 3)+ theme_classic()

########### DENDROGRAM ###################


vsd.distances <- dist(t(assay(vsd.adrenal)))
vsd.hc <- hclust(vsd.distances, method = "ward.D2")
vsd.hc$order
#dend
plot(vsd.hc, hang = -1, cex = 1)
#dend.stage
dendrogram.adrenal.stages <- plot(vsd.hc,labels=stages.adrenal$Stage, hang = -1, cex = 1)

#dend.sample.dev.or
dendrogram.adrenal.label <- plot(vsd.hc,labels=stages.adrenal$sample.dev.origin, hang = -1, cex = 1)


############# SAMPLE DISTANCES ###########################

##### ADRENAL AND CONROLS ###############

sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <-phenodata$Sample_stage
colnames(sampleDistMatrix) <- phenodata$Sample_stage

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.sampleDist.all<- pheatmap(sampleDistMatrix, fontsize = 6,
                                  clustering_distance_rows = sampleDists,
                                  clustering_distance_cols = sampleDists, col=colors)


##### ADRENAL ###############


sampleDists.adrenal <- dist(t(assay(vsd.adrenal)))
sampleDists.adrenal
sampleDistMatrix.adrenal <- as.matrix( sampleDists.adrenal )

rownames(sampleDistMatrix.adrenal) <-stages.adrenal$sample.dev.origin
colnames(sampleDistMatrix.adrenal) <- stages.adrenal$sample.dev.origin
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.tissue.stage.adrenal.labels <- pheatmap(sampleDistMatrix.adrenal, fontsize = 6,
                                         clustering_distance_rows = sampleDists.adrenal,
                                         clustering_distance_cols = sampleDists.adrenal, col=colors)

rownames(sampleDistMatrix.adrenal) <-stages.adrenal$Stage
colnames(sampleDistMatrix.adrenal) <- stages.adrenal$Stage

heatmap.tissue.stage.adrenal <- pheatmap(sampleDistMatrix.adrenal, fontsize = 6,
                                                clustering_distance_rows = sampleDists.adrenal,
                                                clustering_distance_cols = sampleDists.adrenal, col=colors)

dev.off()

#BUILDING RESULTS TABLE. DIFFERENTIAL EXPRESSION ANALYSIS ######### ADRENAL vs CONTROLS ########################

resultsNames(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Plot dispersion estimates
plotDispEsts(dds)

## Define contrasts, extract results table, and shrink the log2 fold changes

contrast <- c("tissue", "Adrenal", "Control")

res <- results(dds, contrast=contrast, alpha = 0.05, independentFiltering=FALSE, cooksCutoff = FALSE)

res_shr <- lfcShrink(dds, contrast=contrast, res=res) # Using type="normal" for shr. Using ashr gives more "strange" results

DESeq2::plotMA(res, ylim=c(-9,9))

DESeq2::plotMA(res_shr, ylim=c(-9,9))



convert <- function(dataset) {
  ensembl = useEnsembl(biomart = "ensembl", host = "www.ensembl.org", mirror='asia')
  mart <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  genes <- rownames(dataset)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  dataset$ensembl <- rownames(dataset)
  dataset <- merge(as.data.frame(dataset),G_list,by.x="ensembl",by.y="ensembl_gene_id")
  return(dataset)
}


res_shr.symbol <- convert(res_shr)
write.csv(res_shr.symbol, file="results/Adrenal-Controls/results.Adrenal.vs.Control.shr.csv",row.names = FALSE)

res_shr.symbol

########## Heatmap #############

padj.cutoff <- 0.05
lfc.cutoff <- 1.5


res_tableOE_tb <- res_shr.symbol%>%
  as_tibble()


normalized_counts.hm <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="ensembl") %>% 
  as_tibble()



## Order results by padj values

top_genes.all <- res_tableOE_tb %>% 
  arrange(padj) %>%
  arrange(desc(log2FoldChange)) %>% drop_na(hgnc_symbol) %>% distinct(hgnc_symbol, .keep_all= TRUE)

top_genes <- left_join(top_genes.all, normalized_counts.hm, by = c('ensembl'))

top_genes <-  top_genes %>% dplyr::select(hgnc_symbol:A835) 
top_genes <- as.data.frame(top_genes)
top_genes <- top_genes[!(top_genes$hgnc_symbol==""), ]

rownames(top_genes) <- top_genes$hgnc_symbol
top_genes <- top_genes %>% dplyr::select(-one_of("hgnc_symbol"))

anno <- phenodata$sample.dev.origin

heatmap.top40.counts.cluster <- pheatmap(top_genes[1:40,], cluster_rows = F, clustering_distance_cols ="euclidean",clustering_method = "complete", cluster_cols = T,scale = "row", fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = T)


heatmap.top500.counts.cluster <- pheatmap(top_genes[1:500,], cluster_rows = F, cluster_cols = T, scale = "row", clustering_distance_cols ="euclidean",clustering_method = "complete",fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = F)


heatmap.top200.counts.cluster <- pheatmap(top_genes[1:200,], cluster_rows = F, cluster_cols = T, scale = "row", clustering_distance_cols ="euclidean",clustering_method = "complete",fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = F)


vsd.df <- as.data.frame(assay(vsd))

vsd.df <- vsd.df %>%
  rownames_to_column(var="ensembl") %>% 
  as_tibble()

top_genes.vsd <- left_join(top_genes.all, vsd.df, by = c('ensembl'))

top_genes.vsd <-  top_genes.vsd %>% dplyr::select(hgnc_symbol:A835) 
top_genes.vsd <- as.data.frame(top_genes.vsd)
top_genes.vsd <- top_genes.vsd[!(top_genes.vsd$hgnc_symbol==""), ]

rownames(top_genes.vsd) <- top_genes.vsd$hgnc_symbol
top_genes.vsd <- top_genes.vsd %>% dplyr::select(-one_of("hgnc_symbol"))

heatmap.top40.vsd.cluster <- pheatmap(top_genes.vsd[1:40,], cluster_rows = F, cluster_cols = T, scale = "row", fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = T)

heatmap.top500.vsd.cluster <- pheatmap(top_genes.vsd[1:500,], cluster_rows = F, cluster_cols = T, scale = "row", fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = F)

heatmap.top200.vsd.cluster <- pheatmap(top_genes.vsd[1:200,], cluster_rows = F, cluster_cols = T, scale = "row", fontsize_row = 6, labels_col = anno, fontsize_col = 6, show_rownames = F)

dev.off()


###################### VOLCANO PLOT ########################

res_tableOE_tb.v <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.0001 & abs(log2FoldChange) >= 2)

volcano.adr.ctrl <- ggplot(res_tableOE_tb.v) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Adrenal vs Control padj<0.0001 & lFC > 2") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,300)) +
  theme(legend.position = "none",plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25))) + theme_classic()


############## REDUCED MODEL LRT #############################

#BUILDING RESULTS TABLE. DIFFERENTIAL EXPRESSION ANALYSIS ######### ADRENAL XY vs XX ########################

resultsNames(dds.adrenal)

res.adrenal.karyotype <- results(dds.adrenal, contrast=c("sex.adrenal", "XY", "XX"), alpha = 0.05, independentFiltering=FALSE, cooksCutoff = FALSE)

res_shr.adrenal.karyotype <- lfcShrink(dds.adrenal, contrast=c("sex.adrenal", "XY", "XX"), res=res.adrenal.karyotype)

DESeq2::plotMA(res.adrenal.karyotype, ylim=c(-9,9))

DESeq2::plotMA(res_shr.adrenal.karyotype, ylim=c(-9,9))

res_shr.adrenal.karyotype.symbol <- convert(res_shr.adrenal.karyotype)

write.csv(res_shr.adrenal.karyotype.symbol, file="results/Adrenal/results.Adrenal.XY.vs.XX.shr.csv", row.names = FALSE)



##############################################

dds.stage <- DESeq(dds.adrenal)

dds_lrt <- DESeq(dds.stage, test="LRT", reduced = ~batch.adrenal + sex.adrenal)

resultsNames(dds_lrt)

normalized_counts.stage <- counts(dds_lrt, normalized=TRUE)
head(normalized_counts.stage)

write.table(normalized_counts.stage, file="results/Adrenal/TimeCourse/normalized_counts.stage.txt", sep="\t", quote=F, col.names=NA)

###### CS23 vs CS20.21 #######

res_LRT_CS23.20 <- results(dds_lrt, contrast=c("stage.adrenal", "CS23", "CS20.21"), independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_CS23.20 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "CS23", "CS20.21"), res=res_LRT_CS23.20)

DESeq2::plotMA(res_LRT_CS23.20, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_CS23.20, ylim=c(-9,9))

res_shr_CS23.20.symbol <- convert(res_LRT_shr_CS23.20)

write.csv(res_shr_CS23.20.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.CS23.vs.CS20.21.shr.csv", row.names = FALSE)

###### F2 vs CS20.21 #######

res_LRT_F2.20 <- results(dds_lrt, contrast=c("stage.adrenal", "F2", "CS20.21"), independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_F2.20 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "F2", "CS20.21"), res=res_LRT_F2.20)

DESeq2::plotMA(res_LRT_F2.20, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_F2.20, ylim=c(-9,9))

res_shr_F2.20.symbol <- convert(res_LRT_shr_F2.20)

write.csv(res_shr_F2.20.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.F2.vs.CS20.21.shr.csv", row.names = FALSE)


###### F4.5 vs CS20.21 #######

res_LRT_F5.20 <- results(dds_lrt, contrast=c("stage.adrenal", "F4.5", "CS20.21"), independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_F5.20 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "F4.5", "CS20.21"), res=res_LRT_F5.20)

DESeq2::plotMA(res_LRT_F5.20, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_F5.20, ylim=c(-9,9))

res_shr_F5.20.symbol <- convert(res_LRT_shr_F5.20)


write.csv(res_shr_F5.20.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.F4.5.vs.CS20.21.csv", row.names = FALSE)

###### F2 vs CS23 #######

res_LRT_F2.23 <- results(dds_lrt, contrast=c("stage.adrenal", "F2", "CS23"),independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_F2.23 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "F2", "CS23"), res=res_LRT_F2.23)

DESeq2::plotMA(res_LRT_F2.23, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_F2.23, ylim=c(-9,9))

res_shr_F2.23.symbol <- convert(res_LRT_shr_F2.23)

write.csv(res_shr_F2.23.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.F2.vs.CS23.shr.csv",row.names = FALSE)

######## F4.5 vs F2 #########

res_LRT_F5.F2 <- results(dds_lrt, contrast=c("stage.adrenal", "F4.5", "F2"), independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_F5.F2 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "F4.5", "F2"), res=res_LRT_F5.F2)

DESeq2::plotMA(res_LRT_F5.F2, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_F5.F2, ylim=c(-9,9))

res_shr_F5.F2.symbol <- convert(res_LRT_shr_F5.F2)

write.csv(res_shr_F5.F2.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.F4.5.vs.F2.shr.csv",row.names = FALSE)


######## F4.5 vs CS23 #########

res_LRT_F5.23 <- results(dds_lrt, contrast=c("stage.adrenal", "F4.5", "CS23"), independentFiltering=FALSE, cooksCutoff = FALSE)

res_LRT_shr_F5.23 <- lfcShrink(dds_lrt, contrast=c("stage.adrenal", "F4.5", "CS23"), res=res_LRT_F5.23)

DESeq2::plotMA(res_LRT_F5.23, ylim=c(-9,9))

DESeq2::plotMA(res_LRT_shr_F5.23, ylim=c(-9,9))

res_shr_F5.23.symbol <- convert(res_LRT_shr_F5.23)

write.csv(res_shr_F5.23.symbol, file="results/Adrenal/TimeCourse/results.Adrenal.F4.5.vs.CS23.shr.csv",row.names = FALSE)

################## Identifying gene clusters exhibiting particular patterns across samples  ################

###### CS23 vs CS20.21 #######

# Subset the LRT results to return genes with padj < 0.05

clustering_sig_genes_CS23.20 <- convert(res_LRT_shr_CS23.20) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.7)

nrow(clustering_sig_genes_CS23.20) #720

clustering_sig_genes_CS23.20 <- clustering_sig_genes_CS23.20[!(clustering_sig_genes_CS23.20$hgnc_symbol==""), ]
nrow(clustering_sig_genes_CS23.20) #607


# Obtain rlog values for those significant genes
vsd.df <- as.data.frame(assay(vsd.adrenal))

cluster_vsd_CS23.20 <- vsd.df[clustering_sig_genes_CS23.20$ensembl, ]

patterns_CS23.20 <- degPatterns(cluster_vsd_CS23.20, metadata = stages.adrenal, time = "Stage")

patterns_CS23.20.genes <- convert(patterns_CS23.20$df)

nrow(patterns_CS23.20.genes) #592


head(res_shr_CS23.20.symbol)
head(patterns_CS23.20.genes)

patterns_CS23.20.genes <- left_join(patterns_CS23.20.genes, res_shr_CS23.20.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)

write.csv(patterns_CS23.20.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.CS23.vs.CS20.21.shr.lFC.0.7.csv", row.names = F)

############ F2 vs CS20.21 ###########


# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F2.20 <- convert(res_LRT_shr_F2.20) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(clustering_sig_genes_F2.20) #1800

clustering_sig_genes_F2.20 <- clustering_sig_genes_F2.20[!(clustering_sig_genes_F2.20$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F2.20) #1613


# Obtain rlog values for those significant genes

cluster_vsd_F2.20 <- vsd.df[clustering_sig_genes_F2.20$ensembl, ]

patterns_F2.20 <- degPatterns(cluster_vsd_F2.20, metadata = stages.adrenal, time = "Stage")

patterns_F2.20.genes <- convert(patterns_F2.20$df)

nrow(patterns_F2.20.genes) #1606

head(res_shr_F2.20.symbol)
head(patterns_F2.20.genes)

patterns_F2.20.genes <- left_join(patterns_F2.20.genes, res_shr_F2.20.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)

write.csv(patterns_F2.20.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.F2.vs.CS20.21.shr.lFC.1.csv", row.names = F)


############ F4.5 vs CS20.21 ###########

# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F5.20 <- convert(res_LRT_shr_F5.20) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(clustering_sig_genes_F5.20) #3603

clustering_sig_genes_F5.20 <- clustering_sig_genes_F5.20[!(clustering_sig_genes_F5.20$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F5.20) #3270


# Obtain rlog values for those significant genes

cluster_vsd_F5.20 <- vsd.df[clustering_sig_genes_F5.20$ensembl, ]

patterns_F5.20 <- degPatterns(cluster_vsd_F5.20, metadata = stages.adrenal, time = "Stage")

patterns_F5.20.genes <- convert(patterns_F5.20$df)

nrow(patterns_F5.20.genes) #3261

head(res_shr_F5.20.symbol)
head(patterns_F5.20.genes)

patterns_F5.20.genes <- left_join(patterns_F5.20.genes, res_shr_F5.20.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)
write.csv(patterns_F5.20.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.F4.5.vs.CS20.21.shr.lFC.1.csv", row.names = F)

###### F2 vs CS23 #######

# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F2.23 <- convert(res_LRT_shr_F2.23) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(clustering_sig_genes_F2.23) #776

clustering_sig_genes_F2.23 <- clustering_sig_genes_F2.23[!(clustering_sig_genes_F2.23$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F2.23) #727


# Obtain rlog values for those significant genes

cluster_vsd_F2.23 <- vsd.df[clustering_sig_genes_F2.23$ensembl, ]

patterns_F2.23 <- degPatterns(cluster_vsd_F2.23, metadata = stages.adrenal, time = "Stage")

patterns_F2.23.genes <- convert(patterns_F2.23$df)

nrow(patterns_F2.23.genes) #719

head(res_shr_F2.23.symbol)
head(patterns_F2.23.genes)

patterns_F2.23.genes <- left_join(patterns_F2.23.genes, res_shr_F2.23.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)
write.csv(patterns_F2.23.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.F2.vs.CS23.shr.lFC.1.csv", row.names = F)

######## F4.5 vs F2 #########

# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F5.F2 <- convert(res_LRT_shr_F5.F2) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(clustering_sig_genes_F5.F2) #323

clustering_sig_genes_F5.F2 <- clustering_sig_genes_F5.F2[!(clustering_sig_genes_F5.F2$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F5.F2) #280


# Obtain rlog values for those significant genes

cluster_vsd_F5.F2 <- vsd.df[clustering_sig_genes_F5.F2$ensembl, ]

patterns_F5.F2 <- degPatterns(cluster_vsd_F5.F2, metadata = stages.adrenal, time = "Stage")

patterns_F5.F2.genes <- convert(patterns_F5.F2$df)

nrow(patterns_F5.F2.genes) #261

head(res_shr_F5.F2.symbol)
head(patterns_F5.F2.genes)

patterns_F5.F2.genes <- left_join(patterns_F5.F2.genes, res_shr_F5.F2.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)
write.csv(patterns_F5.F2.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.F4.5.vs.F2.shr.lFC.1.csv", row.names = F)


######## F4.5 vs CS23 #########

# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F5.23 <- convert(res_LRT_shr_F5.23) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(clustering_sig_genes_F5.23) #2199

clustering_sig_genes_F5.23 <- clustering_sig_genes_F5.23[!(clustering_sig_genes_F5.23$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F5.23) #2027


# Obtain rlog values for those significant genes

cluster_vsd_F5.23 <- vsd.df[clustering_sig_genes_F5.23$ensembl, ]

patterns_F5.23 <- degPatterns(cluster_vsd_F5.23, metadata = stages.adrenal, time = "Stage")

patterns_F5.23.genes <- convert(patterns_F5.23$df)

nrow(patterns_F5.23.genes) #2022

head(res_shr_F5.23.symbol)
head(patterns_F5.23.genes)

patterns_F5.23.genes <- left_join(patterns_F5.23.genes, res_shr_F5.23.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)

write.csv(patterns_F5.23.genes, file="results/Adrenal/TimeCourse/patterns/patterns.Adrenal.F4.5.vs.CS23.shr.lFC.1.csv", row.names = F)


################# SEX DIFFERENCES ##############################
###########################################################################

View(stages.adrenal)
View(countdata.adrenal)

rownames(stages.adrenal) == colnames(countdata.adrenal)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

karyotype.stage.sex <- as.factor(stages.adrenal$stage.sex)

coldata.k <- data.frame(row.names=colnames(countdata.adrenal), batch.adrenal, sex.adrenal, stage.adrenal, karyotype.stage.sex)
coldata.k

dds.k <- DESeqDataSetFromMatrix(countData=countdata.adrenal, colData=coldata.k, design=~batch.adrenal + karyotype.stage.sex)

dds.k <- dds.k[ rowSums(counts(dds.k)) > 1, ]

dds.k <- DESeq(dds.k)
sizeFactors(dds.k)
normalized_counts.k <- counts(dds.k, normalized=TRUE)
head(normalized_counts.k)

resultsNames(dds.k)

res.CS20.21.k <- results(dds.k, contrast=c("karyotype.stage.sex","CS20.21XY","CS20.21XX"), independentFiltering=FALSE, cooksCutoff = FALSE)
res.CS20.21.k.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","CS20.21XY","CS20.21XX"), res = res.CS20.21.k) 
summary(res.CS20.21.k.shr)

res.CS23.k <- results(dds.k, contrast=c("karyotype.stage.sex","CS23XY","CS23XX"), independentFiltering=F, cooksCutoff = FALSE) 
res.CS23.k.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","CS23XY","CS23XX"), res = res.CS23.k) 
summary(res.CS23.k.shr)

res.F2.k <- results(dds.k, contrast=c("karyotype.stage.sex","F2XY","F2XX"), independentFiltering=F, cooksCutoff = FALSE) 
res.F2.k.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","F2XY","F2XX"), res = res.F2.k) 
summary(res.F2.k.shr)

res.F4.5.k <- results(dds.k, contrast=c("karyotype.stage.sex","F4.5XY","F4.5XX"), independentFiltering=F, cooksCutoff = FALSE) 
res.F4.5.k.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","F4.5XY","F4.5XX"), res = res.F4.5.k) 
summary(res.F4.5.k.shr)

res.late.early.XY <- results(dds.k, contrast=c("karyotype.stage.sex","F4.5XY","CS20.21XY"), independentFiltering=F, cooksCutoff = FALSE) 
res.late.early.XY.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","F4.5XY","CS20.21XY"), res = res.late.early.XY) 
summary(res.late.early.XY.shr)

res.late.early.XX <- results(dds.k, contrast=c("karyotype.stage.sex","F4.5XX","CS20.21XX"), independentFiltering=F, cooksCutoff = FALSE) 
res.late.early.XX.shr <- lfcShrink(dds.k, contrast=c("karyotype.stage.sex","F4.5XX","CS20.21XX"), res = res.late.early.XY) 
summary(res.late.early.XX.shr)


res.CS20.21.k.shr.symbol <- convert(res.CS20.21.k.shr)
res.CS23.sex.k.symbol <- convert(res.CS23.k.shr)
res.F2.sex.k.symbol <- convert(res.F2.k.shr)
res.F4.5.sex.k.symbol <- convert(res.F4.5.k.shr)
res.late.early.XY.symbol <- convert(res.late.early.XY.shr)
res.late.early.XX.symbol <- convert(res.late.early.XX.shr)


write.csv(res.CS20.21.k.shr.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.karyotype.shr.CS20.21.XY.vs.XX.csv",row.names = F)
write.csv(res.CS23.sex.k.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.karyotype.shr.CS23.XY.vs.XX.csv",row.names = F)
write.csv(res.F2.sex.k.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.karyotype.shr.F2.XY.vs.XX.csv",row.names = F)
write.csv(res.F4.5.sex.k.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.karyotype.shr.F4.5.XY.vs.XX.csv",row.names = F)
write.csv(res.late.early.XY.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.shr.F4.5.XY.vs.CS20.21.XY.csv",row.names = F)
write.csv(res.late.early.XX.symbol, file="results/Adrenal/TimeCourse/SexDifferences/results.Adrenal.shr.F4.5.XX.vs.CS.20.21.XX.csv",row.names = F)



############ F4.5 XY vs CS20.21 XY ###########


# Subset the LRT results to return genes with padj < 0.05
clustering_sig_genes_F5.20.XY <- convert(res.late.early.XY.shr) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.7)

nrow(clustering_sig_genes_F5.20.XY) #5712

clustering_sig_genes_F5.20.XY <- clustering_sig_genes_F5.20.XY[!(clustering_sig_genes_F5.20.XY$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F5.20.XY) #5340


# Obtain values for those significant genes

cluster_vsd_F5.20.XY <- vsd.df[clustering_sig_genes_F5.20.XY$ensembl, ]

patterns_F5.20.XY <- degPatterns(cluster_vsd_F5.20.XY, metadata = stages.adrenal, time = "Stage")

patterns_F5.20.genes.XY <- convert(patterns_F5.20.XY$df)

nrow(patterns_F5.20.genes.XY) #5331

patterns_F5.20.genes.XY <- left_join(patterns_F5.20.genes.XY, res.late.early.XY.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)

write.csv(patterns_F5.20.genes.XY, file="results/Adrenal/TimeCourse/SexDifferences/TimeCourse/patterns.Adrenal.F4.5.XY.vs.CS20.21.XY.LFC.0.7.shr.csv")


############ F4.5 XX vs CS20.21 XX ###########


# Subset the LRT results to return genes with padj < 0.05

clustering_sig_genes_F5.20.XX <- convert(res.late.early.XX.shr) %>%
  as_tibble() %>% 
  arrange(padj) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.7)

nrow(clustering_sig_genes_F5.20.XX) #3449

clustering_sig_genes_F5.20.XX <- clustering_sig_genes_F5.20.XX[!(clustering_sig_genes_F5.20.XX$hgnc_symbol==""), ]
nrow(clustering_sig_genes_F5.20.XX) #3246


# Obtain values for those significant genes

cluster_vsd_F5.20.XX <- vsd.df[clustering_sig_genes_F5.20.XX$ensembl, ]

patterns_F5.20.XX <- degPatterns(cluster_vsd_F5.20.XX, metadata = stages.adrenal, time = "Stage")

patterns_F5.20.genes.XX <- convert(patterns_F5.20.XX$df)

nrow(patterns_F5.20.genes.XX) #3229

patterns_F5.20.genes.XX <- left_join(patterns_F5.20.genes.XX, res.late.early.XX.symbol, by = c('ensembl','hgnc_symbol')) %>% dplyr::select(-genes)

write.csv(patterns_F5.20.genes.XX, file="results/Adrenal/TimeCourse/SexDifferences/TimeCourse/patterns.Adrenal.F4.5.XX.vs.CS20.21.XX.LFC.0.7.shr.csv")


################ TPM PLOTS #################

######### ADRENAL & CONTROLS ############

head(normalized_counts)

gene.length <- read.table("data/genelength.csv", row.names=1, sep=",",header=T)
head(gene.length)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

dim(normalized_counts)
rownames(normalized_counts)
rownames(gene.length)
dim(gene.length)

gene.length.all <- subset(gene.length, (rownames(gene.length) %in% rownames(as.data.frame(normalized_counts))))

tpms.all <- apply(normalized_counts, 2, function(x) tpm(x, gene.length.all$Length))

colSums(tpms.all)
colMeans(tpms.all)

write.table(tpms.all, file="results/TPM.All.csv", sep=",")

head(tpms.all)


######## ADRENAL ############

head(normalized_counts.adrenal)

gene.length <- subset(gene.length, (rownames(gene.length) %in% rownames(normalized_counts.adrenal)))

tpms <- apply(normalized_counts.adrenal, 2, function(x) tpm(x, gene.length$Length))

colSums(tpms)
colMeans(tpms)

write.table(tpms, file="results/TPM.Adrenal.csv", sep=",")

colnames(tpms) == rownames(stages.adrenal)


########## Heatmap ADRENAL STEROIDOGENESIS #############

head(normalized_counts.adrenal)
head(vsd.df)

counts_heatmap <- as.data.frame(normalized_counts.adrenal)
head(counts_heatmap)

counts_heatmap <- convert(counts_heatmap)
nrow(counts_heatmap)

counts_heatmap <- counts_heatmap[!(counts_heatmap$hgnc_symbol==""), ]

counts_heatmap <- counts_heatmap[!duplicated(counts_heatmap$hgnc_symbol),]

rownames(counts_heatmap) <- counts_heatmap$hgnc_symbol

counts_heatmap <- counts_heatmap %>% dplyr::select(-one_of("ensembl","hgnc_symbol"))
dim(counts_heatmap)

steroido <- c("MC2R", "MRAP", "STAR", "CYP11A1", "HSD3B2", "CYP17A1", "POR", "CYB5A", "SULT2A1", "PAPSS2","CYP21A2", "CYP11B1", "CYP11B2")
tfs <- c("FOSL1", "MAFF", "ZNF774", "ARX", "HOXA5", "NR4A2", "NR4A3", "HES6", "KLF4", "ZNF331", "NR4A1", "SIX2", "BHLHE41")

hmap <- counts_heatmap[steroido,]
hmap.tf <- counts_heatmap[tfs,]
hmap

heatmap.steroido.counts.cluster <- pheatmap(hmap, cluster_rows = F, cluster_cols = T, scale = "row", labels_col = sample.adrenal.dev.origin)

hmap.tf.counts.cluster <- pheatmap(hmap.tf, cluster_rows = F, cluster_cols = T, scale = "row",labels_col = sample.adrenal.dev.origin)

dev.off()

head(vsd.df)

counts_heatmap.vsd <- vsd.df
head(counts_heatmap.vsd)

counts_heatmap.vsd <- convert(counts_heatmap.vsd)
nrow(counts_heatmap.vsd)

counts_heatmap.vsd <- counts_heatmap.vsd[!(counts_heatmap.vsd$hgnc_symbol==""), ]

counts_heatmap.vsd <- counts_heatmap.vsd[!duplicated(counts_heatmap.vsd$hgnc_symbol),]

rownames(counts_heatmap.vsd) <- counts_heatmap.vsd$hgnc_symbol

counts_heatmap.vsd <- counts_heatmap.vsd %>% dplyr::select(-one_of("ensembl","hgnc_symbol"))
dim(counts_heatmap.vsd)


steroido <- c("MC2R", "MRAP", "STAR", "CYP11A1", "HSD3B2", "CYP17A1", "POR", "CYB5A", "SULT2A1", "CYP21A2", "CYP11B1", "CYP11B2")
hmap.vsd <- counts_heatmap.vsd[steroido,]
hmap.vsd.tf <- counts_heatmap.vsd[tfs,]

hmap.vsd

heatmap.steroido.vsd.cluster <- pheatmap(hmap.vsd, cluster_rows = F, cluster_cols = T, scale='row', labels_col = sample.adrenal.dev.origin)
heatmap.vsd.tf.cluster <- pheatmap(hmap.vsd.tf, cluster_rows = F, cluster_cols = T, scale = "row", labels_col = sample.adrenal.dev.origin)

dev.off()



########################################## FUNCTIONS #####################################
##########################################           #####################################

########### Counts Plot ##############

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Obtained from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


save.image(file = "Adrenal.RData")

#############################################################################
#############################################################################





colData(dds) # dds- all samples
colData(dds.adrenal) # dds.adrenal Only adrenal samples



counts_plot <- function(g){
  geneName <- mapIds(Homo.sapiens,keys=g,column="ENSEMBL",keytype="SYMBOL")
  geneCounts <- plotCounts(dds, gene = geneName, intgroup = c("batch", "tissue", "sex", "stage", "origin"),returnData = TRUE)
  t <- ggplot(geneCounts, aes(x = tissue, y = count)) + geom_violin(trim = T,aes(fill=tissue))  + theme_classic() + scale_y_log10() + ggtitle(g)+ geom_boxplot(fill = "grey", width=0.1)+ scale_fill_brewer(palette="Dark2")
  k <- ggplot(geneCounts, aes(x = tissue, y = count)) + geom_violin(trim = T, aes(fill=tissue))  + theme_classic() + scale_y_log10() + ggtitle(g) + facet_wrap(~sex)+ geom_boxplot(fill = "grey", width=0.1)+ geom_boxplot(fill = "grey", width=0.1)+ scale_fill_brewer(palette="Dark2")
  s <- ggplot(geneCounts, aes(x = tissue, y = count)) + geom_violin(trim = T, aes(fill=tissue))  + theme_classic() + scale_y_log10() + ggtitle(g) + facet_wrap(~stage)+ geom_boxplot(fill = "grey", width=0.1)+ scale_fill_brewer(palette="Dark2")
  print(geneCounts)
  #list(t,k,s) # If we need graph on same figure use multiplot(t,k, cols = 2)
  multiplot(t,k,s, cols = 2)
}


counts_plot("ACE2")
counts_plot("CCN3")

counts_plot("NEUROD4")
counts_plot("FABP1")
counts_plot("INA")
counts_plot("APOA2")

counts_plot("NR5A1")
counts_plot("IL1RL1")


colData(dds.adrenal)

counts_plot_sex <- function(g){
  geneName <- mapIds(Homo.sapiens,keys=g,column="ENSEMBL",keytype="SYMBOL")
  geneCounts <- plotCounts(dds.adrenal, gene = geneName, intgroup = c("batch.adrenal", "stage.sex"),returnData = TRUE)
  t <- ggplot(geneCounts, aes(x = stage.sex, y = count, fill = stage.sex)) + theme_classic() + scale_y_log10() + ggtitle(g)+ geom_violin(trim = T)+ geom_boxplot(fill = "grey", width=0.1)
  print(geneCounts)
  list(t) # If we need graph on same figure use multiplot(t,k, cols = 2)
}


counts_plot_sex("ACE2")
counts_plot_sex("CCN3")


counts_plot.adrenal <- function(g){
  geneName <- mapIds(Homo.sapiens,keys=g,column="ENSEMBL",keytype="SYMBOL")
  geneCounts <- plotCounts(dds.adrenal, gene = geneName, intgroup = c("batch.adrenal", "sex.adrenal", "stage.sex", "stage.adrenal"),returnData = TRUE)
  t <- ggplot(geneCounts, aes(x = stage.adrenal, y = count)) + geom_violin(trim=T,aes(fill=stage.adrenal))  + theme_classic() + scale_y_log10() + geom_boxplot(fill = "grey",width=0.1) +
    ggtitle(g)
  p <- ggplot(geneCounts, aes(x = stage.adrenal, y = count)) + geom_violin(trim=T,aes(fill=stage.adrenal))  + theme_classic() + scale_y_log10() + geom_boxplot(fill = "grey", width=0.1) +
    ggtitle(g) + facet_wrap(~sex.adrenal)
  print(geneCounts)
  #print(t)
  #list(t,p) # If we need graph on same figure use multiplot(t,k, cols = 2)
  multiplot(t,p, cols = 2)
}


counts_plot.adrenal("CD99L2")

counts_plot.adrenal("ACE2")


tpms_plot <- function(g){
  geneName <- mapIds(Homo.sapiens,keys=g,column="ENSEMBL",keytype="SYMBOL")
  gene <- as.data.frame(tpms[geneName,])
  gene$Stage <- stages.adrenal$Stage
  gene$Sex <- stages.adrenal$Sex
  gene$Stage.Sex <- stages.adrenal$stage.sex
  print(gene)
  p <- ggplot(gene,aes(x=Stage ,y = tpms[geneName, ], fill=Stage))+geom_violin(trim = T)+ geom_boxplot(fill = "grey", width=0.1)+facet_wrap(~Sex) + labs (y = "TPM") + ggtitle(g) +theme_classic()
  q <- ggplot(gene,aes(x=Stage.Sex ,y = tpms[geneName, ], fill=Stage))+ labs (y = "TPM") + ggtitle(g)+theme_classic()+geom_violin(trim = T)+ geom_boxplot(fill = "grey", width=0.1)
  r <- ggplot(gene,aes(x=Stage ,y = tpms[geneName, ], fill=Stage))+ labs (y = "TPM") + ggtitle(g)+theme_classic()+geom_violin(trim = T)+ geom_boxplot(fill = "grey", width=0.1)
  list(p,q,r)
}

tpms_plot("CD99L2")



dev.off()








