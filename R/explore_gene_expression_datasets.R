#Script to obtain microarray and RNAseq data listed in the ALS PandaOmics
#paper from GEO and ArrayExpress

#Libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(GEOquery)
library(ggplot2)
library(DESeq2)
options(stringsAsFactors = F)

#Create a vector with all GEO Ids for the cited datasets
table_01 <- fread("data/table_01.txt")
gse_ids <- table_01 %>%
  filter(grepl("GSE",Data_series)) %>%
  pull(Data_series) %>%
  unique()

####GSE67196 - fALS and sALS - Frontal Cortex and Cerebellum####
#Download GSE data
gse <- getGEO(gse_ids[1])

#Download raw count data from supplementary files
count_data <- fread(gse$GSE67196_series_matrix.txt.gz@experimentData@other$supplementary_file)
count_data <- count_data %>%
  dplyr::select(-1)

#Obtain metadata
GEO_metadata <- as.data.frame(gse$GSE67196_series_matrix.txt.gz) %>%
  mutate(title=gsub(x = title,pattern = "cereb_a|cereb_b",replacement = "cereb"))


#Harmonize metadata sample ids with count data colnames
samples <- data.frame(sample=colnames(count_data)[-1]) %>%
  mutate(title=gsub(sample,pattern = "ALS0|ALS00|ALS",replacement = ""),
         title=gsub(title,pattern = "fcx",replacement = "FCX"))

#Create pheno data
pheno_data <- GEO_metadata %>%
  dplyr::left_join(samples,by = "title") %>%
  dplyr::select(41,39,40) %>%
  dplyr::rename(genotype=2,brain_region=3) %>%
  dplyr::mutate(brain_region=gsub(x = brain_region,pattern = " ",replacement = "_"),
                genotype=factor(x = genotype,levels = c("Healthy",
                                                        "c9ALS",
                                                        "sALS")),
                group=paste0(brain_region,"_",genotype))

#Run DESeq2 for each comparisson
#Remove duplicated gene names and Excell errors
count_data_clear <- count_data %>%
  filter(!grepl("-Mar|-Sep|-Dec",GeneID),
         !duplicated(GeneID)) %>%
  column_to_rownames("GeneID")

count_data_in <- count_data_clear[,pheno_data$sample]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data_in,
                                      colData = pheno_data,
                                      design = ~group,
                                      tidy = F)

#Remove lowly expressed genes
keep_10 <- rowSums(DESeq2::counts(dds) >= 10) >= 10
dds <- dds[keep_10,]

#PCA of samples
vsd <- DESeq2::vst(dds,blind = F)

#Calculate and plot PCA manually (same way DESeq2 does it)
# calculate the variance for each gene
rv <- rowVars(assay(vsd))
# select the ntop genes by variance
ntop=500
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

vsd_500 <- assay(vsd)[select,]
pca <- prcomp(t(vsd_500))

pc1_variance <- round(100*(pca$sdev[1]^2/sum(pca$sdev^2)))
pc2_variance <- round(100*(pca$sdev[2]^2/sum(pca$sdev^2)))
pca_df_plot <- pca$x %>%
  as.data.frame() %>%
  select(c(1,2)) %>%
  rownames_to_column("sample") %>%
  left_join(pheno_data,by = "sample") %>%
  column_to_rownames("sample")

p <- pca_df_plot %>%
  mutate(group=paste0(brain_region,":",genotype)) %>%
  ggplot(aes(x=PC1,y = PC2,color=genotype,pch=brain_region))+
  geom_point(size=5)+
  scale_color_brewer(palette = "Set1")+
  scale_x_continuous(name = paste0("PC1: ",pc1_variance,"% variance"))+
  scale_y_continuous(name = paste0("PC2: ",pc2_variance,"% variance"))+
  theme_bw()
p

#Run DE analysis before batch removal
dds <- DESeq2::DESeq(object = dds)

#Remove hidden batch effect(s)
#(PCA shows there is an unknown variable that splits the data with
#5% of variance explanation and also a very large difference
#between clusters of samples in PC1)
library(sva)
dds <- estimateSizeFactors(dds)
dat <- DESeq2::counts(object = dds,normalized=T)
idx  <- rowSums(DESeq2::counts(dds) >= 10) >= 10
dat  <- dat[idx, ]
mod  <- model.matrix(~group,colData(dds))
mod0 <- model.matrix(~1, colData(dds))
n.sv <- sva::num.sv(dat = dat,mod = mod)
svseq <- svaseq(dat, mod, mod0,n.sv = 2)

sva_results <- svseq$sv

#Plot new PCA with removed surrogate variables
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleanMX <- cleanY(y = dat,mod = mod,svs = sva_results)
# calculate the variance for each gene
rv <- rowVars(cleanMX)
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

pca <- prcomp(t(cleanMX[select,]))

pc1_variance <- round(100*(pca$sdev[1]^2/sum(pca$sdev^2)))
pc2_variance <- round(100*(pca$sdev[2]^2/sum(pca$sdev^2)))
pca_df <- pca$x

pca_df_plot <- pca_df %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(pheno_data,by = "sample")

p <- pca_df_plot %>%
  mutate(group=paste0(brain_region,":",genotype)) %>%
  ggplot(aes(x=PC1,y = PC2,color=group,pch=brain_region))+
  geom_point(size=5)+
  scale_color_brewer(palette = "Set1")+
  scale_x_continuous(name = paste0("PC1: ",pc1_variance,"% variance"))+
  scale_y_continuous(name = paste0("PC2: ",pc2_variance,"% variance"))+
  theme_bw()
p

#Run DESeq2 removing 2 surrogate variables
#(and not removing for comparison)
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# ddssva$SV4 <- svseq$sv[,4]
design(ddssva) <- ~ SV1 + SV2 + group

ddssva <- DESeq2::DESeq(object = ddssva)

#Get results and plot volcano plots
conditions <- c("c9ALS","sALS")
regions <- c("Cerebellum","Frontal_Cortex")
contrast_df <- expand.grid(conditions,regions) %>%
  dplyr::rename(condition=1,region=2) %>%
  dplyr::mutate(A=paste0(region,"_",condition),
                B=paste0(region,"_","Healthy"))

for(i in 1:nrow(contrast_df)){
  contr <- c("group",contrast_df$A[i],contrast_df$B[i])
  res_sva <- as.data.frame(results(ddssva,contrast = contr)) %>%
    mutate(model="sva",
           region=contrast_df$region[i],
           comparison=paste0(contrast_df$condition[i],"_vs_Healthy")) %>%
    rownames_to_column("Gene")
  
  res_dds <- as.data.frame(results(dds,contrast = contr)) %>%
    mutate(model="original",
           region=contrast_df$region[i],
           comparison=paste0(contrast_df$condition[i],"_vs_Healthy")) %>%
    rownames_to_column("Gene")
  
  res_out <- rbind(res_dds,res_sva)
  if(i==1){
    res_all <- res_out
  }else{
    res_all <- rbind(res_all,res_out)
  }
}

plot_df <- res_all %>%
  filter(!is.na(padj)) %>%
  mutate(significativo=ifelse(pvalue<0.05 & abs(log2FoldChange) >= 1,
                              yes = "yes",
                              no = "no"),
         dir=ifelse(log2FoldChange>0,yes = "positive",no = "negative"),
         DE=ifelse(significativo=="yes" & dir=="positive",
                   yes = "positive",
                   no = ifelse(significativo=="yes" & dir=="negative",
                               yes = "negative",
                               no = "not signficant")),
         DE=factor(DE, levels = c("negative","not signficant","positive")))

ndegs <- plot_df %>%
  filter(significativo=="yes") %>%
  mutate(number=1) %>%
  group_by(region,comparison,model,dir) %>%
  summarise(nDEGs=sum(number)) %>%
  ungroup()


plot <- plot_df %>%
  ggplot(aes(x=log2FoldChange,y = -log(pvalue))) +
  geom_point(aes(color=DE)) +
  scale_color_manual(values = c("blue","grey","red3")) +
  geom_hline(yintercept = -log(0.05),linetype="dashed") +
  geom_vline(xintercept = c(-1,1),linetype="dashed") +
  theme_bw()+
  facet_wrap(facets = ~region+comparison+model,
             ncol = 4,nrow = 2)
plot
