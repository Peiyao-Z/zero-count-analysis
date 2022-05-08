rm(list=ls())
library(SingleCellExperiment)
library(pbmcapply)
library(tidyverse)
library(fields)
library(Seurat)
library(CARD)
library(parallel)
library(foreach)

datapath = "FILE_DATAPATH"

# ST-HBC ---------
spatial_count <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
spatial_count = t(spatial_count)

CARDfree_obj = readRDS('ST_HBC_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'ST_MOB_ref_sum.rds')
# ST-HBC ---------
spatial_count <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
spatial_count = t(spatial_count)

CARDfree_obj = readRDS('ST_HBC_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'ST_HBC_ref_sum.rds')


## 10x_MB(S-A)----------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrainAnterior10x.h5"))

CARDfree_obj = readRDS('10x_MB(S-A)_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_MB(S-A)_ref_sum.rds')

# 10x_MB(S-P)----------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrainPosterior10x.h5"))

CARDfree_obj = readRDS('10x_MB(S-P)_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_MB(S-P)_ref_sum.rds')

## 10x_MB(C)----------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrain10x.h5"))

CARDfree_obj = readRDS('10x_MB(C)_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_MB(C)_ref_sum.rds')

## 10x_HH----------
spatial_count<-Read10X_h5(paste0(datapath,"./HumanHeart10x.h5"))

CARDfree_obj = readRDS('10x_HH_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_HH_ref_sum.rds')

## 10x_HBC----------
spatial_count<-Read10X_h5(paste0(datapath,"./BC10x.h5"))

CARDfree_obj = readRDS('10x_HBC_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_HBC_ref_sum.rds')


## 10x_MK----------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseKidney10x.h5"))

CARDfree_obj = readRDS('10x_MK_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'10x_MK_ref_sum.rds')

## 10x_Slideseq----------
spatial_count<-readRDS(paste0(datapath,"Slide-seq.rds"))

CARDfree_obj = readRDS('Slideseq_ref.rds')
cov = CARDfree_obj@Proportion_CARD
tmp = matrix(0,nrow = length(setdiff(colnames(spatial_count),rownames(cov))),ncol = ncol(cov))
rownames(tmp) = setdiff(colnames(spatial_count),rownames(cov))
cov = rbind(cov, tmp)
cov = cov[match(colnames(spatial_count),rownames(cov)),]

source("CARD_overdispersion.R")

saveRDS(sum_stat,'Slideseq_ref_sum.rds')



