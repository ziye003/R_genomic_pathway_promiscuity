# Databricks notebook source
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('scater')
library(scater)

# COMMAND ----------

install.packages(c('R.utils','mclust'))

# COMMAND ----------

library(data.table)
library(corrr)
library(mclust)
library(MASS)
library(scater)

# COMMAND ----------

# MAGIC %md # load unfiltered data

# COMMAND ----------

df=fread('/dbfs/mnt/client-002sap21p010-pasca/04_data_analysis/temp/GSM3015987_H20180123_neuron.csv')
df[1:3,1:4]

# COMMAND ----------

Dates=fread('/dbfs/mnt/client-002sap21p010-pasca/04_data_analysis/temp/dates_neuron.csv')
Dates[1:4]
Date=as.data.frame(unlist(Dates))
unique(Date)

# COMMAND ----------

Stage=Date
Stage[which(Date=="naiveH1"),]=1
Stage[which(Date=="primeH1"),]=2
Stage[which(Date=="primeH9"),]=2
Stage[which(Date=="eb1"),]=3
Stage[which(Date=="eb2"),]=4
Stage[which(Date=="eb3"),]=5
Stage[which(Date=="eb4"),]=6
Stage[which(Date=="ebD4"),]=7

# COMMAND ----------

# MAGIC %md # preprocess

# COMMAND ----------

df=as.data.frame(df)
rownames(df)=df$V1
df = df[,-c(1)]

# COMMAND ----------

impute_missing <- function(x) {
  miss_filter <- is.na(x) | is.nan(x) | is.null(x) | x == 0
  if (sum(miss_filter) == 0) return(x)
  min_val <- min(x[!miss_filter])
  x[miss_filter] <- stats::runif(sum(miss_filter), min = min_val/10, max = min_val)
  return(x)
}

scale_conti = function(x){
  sdx = sd(x,na.rm=TRUE)
  newx = x/sd(x)
  return(newx)
}

# COMMAND ----------

# select the top 1000 most variant gene the gene by variance 
df$variance <- apply(df, 1, var)
df=df[order(df$variance, decreasing = T),]
df = subset(df, select = -variance )
print(dim(df))
top_1000_df<-df[1:1000,]

pos_scale_df=sapply(top_1000_df, function(x) scale_conti(x))
pos_scale_df=as.data.frame(pos_scale_df)
                    
pos_scale_df[1:3,1:4]
print(dim(pos_scale_df))

# COMMAND ----------

# MAGIC %md # nmf cluster

# COMMAND ----------

# Install
install.packages('NMF')
# Load
library(NMF)

# COMMAND ----------

# find best number of cluster 
estim.r <- nmf(pos_scale_df,4:7, nrun=2, seed=123456)
plot(estim.r)

# COMMAND ----------

res <- nmf(pos_scale_df,5,nrun=2)
res
ClusterID=(as.matrix(apply(t(coef(res)),1,which.max)))
ClusterID
dim(ClusterID)

saveRDS(ClusterID, '/dbfs/mnt/client-002sap21p010-pasca/04_data_analysis/temp/hN5_NMF.RData')
ClusterID=readRDS('/dbfs/mnt/client-002sap21p010-pasca/04_data_analysis/temp/hN5_NMF.RData')

# COMMAND ----------

Cluster7ID=readRDS('/dbfs/mnt/client-002sap21p010-pasca/04_data_analysis/temp/hN7_NMF.RData')

# COMMAND ----------

# MAGIC %md # pca

# COMMAND ----------

# install.packages('ggfortify')
library(ggfortify)

# COMMAND ----------

df <- t(pos_scale_df)
pca_res <- prcomp(df, scale. = TRUE)

# COMMAND ----------

Date=as.data.frame(unlist(Date))
Stage=as.data.frame(unlist(Stage))
Cluster7ID=as.data.frame(as.character(unlist(Cluster7ID)))
RN=as.data.frame(as.character(unlist(colnames(pos_scale_df))))

colnames(Date)='Date'
colnames(Stage)='Stage'
colnames(Cluster7ID)='Cluster7ID'
colnames(RN)='RN'

# COMMAND ----------

all=cbind(df,Stage,Date,Cluster7ID,RN)
autoplot(pca_res, data=all,colour = 'Cluster7ID')

# COMMAND ----------

autoplot(pca_res, data=all,colour = 'Date')

# COMMAND ----------

# MAGIC %md # tsne

# COMMAND ----------

library(tidyverse)
library(Rtsne)
theme_set(theme_bw(18))

# COMMAND ----------

tsne <- all %>% 
  drop_na() %>%
  mutate(ID=row_number()) 


# COMMAND ----------

tsne_meta <- all %>%
  select(Date,Stage,Cluster7ID,RN)

# COMMAND ----------

set.seed(142)
rownames(tsne) <- NULL
tSNE_fit <- tsne %>%
  select(where(is.numeric)) %>%
  column_to_rownames(RN) %>%
  scale() %>% 
  Rtsne()

# tSNE_fit[1:3,1:4]
print(dim(tSNE_fit$Y))

# COMMAND ----------

tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())
rownames(tSNE_df)=unlist(RN)
tSNE_df['RN']=unlist(RN)
tSNE_df[1:4,1:2]

# COMMAND ----------

tSNE_df <- tSNE_df %>%
  inner_join(tsne_meta, by="RN")

# COMMAND ----------

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Cluster7ID))+
  geom_point()+
  theme(legend.position="bottom")
# ggsave("tSNE_plot_example1.png")

# COMMAND ----------

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Date))+
  geom_point()+
  theme(legend.position="bottom")
# ggsave("tSNE_plot_example1.png")

# COMMAND ----------

# MAGIC %md #umap

# COMMAND ----------

# install.packages('umap')

# COMMAND ----------

umap_fit <- tsne %>%
  select(where(is.numeric)) %>%
#   column_to_rownames(RN) %>%
  scale() %>% 
  umap(n_components = 2, random_state = 15)
umap_fit

# COMMAND ----------

umap_df <- umap_fit$layout %>% 
  as.data.frame() %>%
  rename(umap1="V1",
         umap2="V2") %>%
  mutate(ID=row_number())
rownames(umap_df)=unlist(RN)
umap_df['RN']=unlist(RN)
umap_df[1:4,1:2]

# COMMAND ----------

umap_df <- umap_df %>%
  inner_join(tsne_meta, by="RN")

# COMMAND ----------

umap_df %>%
  ggplot(aes(x = umap1, 
             y = umap2,
             color = Cluster7ID))+
  geom_point()+
  theme(legend.position="bottom")

# COMMAND ----------

umap_df %>%
  ggplot(aes(x = umap1, 
             y = umap2,
             color = Date))+
  geom_point()+
  theme(legend.position="bottom")

# COMMAND ----------


