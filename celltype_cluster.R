library(parallel)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(MASS)
require(pscl)
library(foreach)
library(doParallel)
library(tidyverse)
library(vcdExtra)
library(lmtest)
library(Seurat)
library(Matrix)

datapath = "datapath"

clusterPC<-function(tmp,npcs = 15,resolution =0.5){
  colnames(data)<-c(1:ncol(data))
  pbmc <- CreateSeuratObject(counts = data)
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  levels(Idents(pbmc))
  return(pbmc)
}

source("function_offset.R")

## Slide-seq2 ----
data<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
tmp<-clusterPC(tmp = data)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/Slide-seqV2_celltype.rds")
rm(data,tmp,cellType_sum)


# ## Slide-seq -----
data<-readRDS(paste0(datapath,"Slide-seq.rds"))
#cluster<-read.csv("D:/droplet/Slideseq/AnalogizerClusterAssignments.csv",header = T)
cluster<-read.csv(paste0(datapath,"AnalogizerClusterAssignments.csv"),header = T)
#cluster<-cluster[,c("bc","cluster")]
group<-unique(cluster$x)
cellType_sum <- list()
for (cellType in 1:11) {
  tt<-cluster$Var1[cluster$x==cellType]
  dat<-data[,colnames(data)%in%tt]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  print(group)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info,tt)
}
saveRDS(cellType_sum,"celltype/Slide-seq_celltype.rds")
rm(data,dat,tmp,cellType_sum)


# Paired_cell_sequencing ----
data<-read.csv(paste0(datapath,"liver single cell.txt"),sep = "\t")
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/pcRNAseq_celltype.rds")
rm(data,dat,tmp,cellType_sum)


## LCM -----
data<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
cluster<-read.csv(paste0(datapath,"cell_zone_table.txt"),sep = "")
colnames(data)<-cluster$zone
#group<-levels(as.factor(cluster$zone))
group<-names(sort(table(cluster$zone),decreasing = T))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
  print(group)
}
saveRDS(cellType_sum,"celltype/LCM_celltype.rds")
rm(data,dat,tmp,cellType_sum)


# ##seqFISH -----
library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==43)]
data = dat
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("D:/droplet/gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
  print(group)
}
saveRDS(cellType_sum,"celltype/seqFISH_celltype.rds")
rm(data,dat,tmp,cellType_sum)

## HDST ----
load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
data<-sp_count
rm(location,sp_count)
ann = read.csv(paste0(datapath,"annotation.csv"))
ann$new = paste0(ann$spot_x,"x",ann$spot_y)
# ttt = read.table("CN24_D1_HDST_cell_type_assignments_final.tsv")
# cluster<-read.table(paste0(datapath,"CN24_D1_segments_cell_type_assignments.tsv"))
# group<-unique(cluster$cell_type)
group = unique(ann$poly.ID)
cellType_sum <- list()
for (cellType in group) {
  tt = ann$new[which(ann$poly.ID==cellType)]
  #tt<-cluster$bc[which(cluster$cell_type==cellType)]
  dat <- data[,colnames(data)%in%tt]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  print(group)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info,tt)
}
saveRDS(cellType_sum,"celltype/HDST_celltype.rds")
rm(data,tmp,cellType_sum)

data<-read.csv(paste0(datapath,"MERFISH.csv"))
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
data<-t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/MERFISH_celltype.rds")
rm(data,dat,tmp,cellType_sum)

# STARmap----
data<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
data<-t(data)
cluster<-read.csv(paste0(datapath,"class_labels.csv"))
#cluster<-read.csv("D:/droplet/STARmap/class_labels.csv")
clu<-data.frame(CellID = 0:(ncol(data)-1))
cluster<-left_join(clu,cluster,by="CellID")
colnames(data)<-cluster$ClusterName
group<-names(sort(table(cluster$ClusterName),decreasing = T))

cellType_sum <- list()

for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/STARmap_celltype.rds")
rm(data,dat,tmp,cellType_sum)

## ST----
data <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/ST_BC_celltype.rds")
rm(data,dat,tmp,cellType_sum)

data <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/ST_MOB_celltype.rds")
rm(data,dat,tmp,cellType_sum)


## 10x-----------
#AdultMouseBrain10x
data<-readRDS(paste0(datapath,"MouseBrain10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_MouseBrain(Coronal)_celltype.rds")
rm(data,dat,tmp,cellType_sum)

#BC
data<-readRDS(paste0(datapath,"BC10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_HBC_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-readRDS(paste0(datapath,"HumanHeart10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_HH_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-readRDS(paste0(datapath,"HumanLymph10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_HL_celltype.rds")
rm(data,dat,tmp,cellType_sum)



data<-readRDS(paste0(datapath,"MouseKidney10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_MK_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-readRDS(paste0(datapath,"MouseBrainAnterior10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_MB(S-A)_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-readRDS(paste0(datapath,"MouseBrainPosterior10x_count.rds"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/10x_MB(S-P)_celltype.rds")
rm(data,dat,tmp,cellType_sum)



## tomo -----
data<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(tmp = data)
colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype/Tomo-seq_celltype.rds")
rm(data,dat,tmp,cellType_sum)


## seqFISH -----
# #setwd("D:/droplet/seqFISH")
data<-read.csv(paste0(datapath,"ob_counts.csv"))
data<-t(data)
cluster<-read.csv(paste0(datapath,"OB_cell_type_annotations.csv"))
cellType_sum <- list()
for (cellType in 1:25) {
  #print(cellType)
  dat<-data[,which(cluster$louvain==cellType)]
  #dat <- data[,colnames(data)%in%cluster$louvain]
  source("gene-wise_offset.R")
  cellType_sum[[cellType]] <- sum_stat
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype new/seqFISH+_celltype.rds")
rm(data,dat,tmp,cellType_sum)




