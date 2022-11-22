rm(list=ls())
library(SingleCellExperiment)
library(pbmcapply)
library(CARD)
library(tidyverse)
library(fields)
library(Seurat)
library(tidyr)
library(dplyr)
datapath = "ST_file_DATAPATH"
DATAPATH = "scRNA_file_DATAPATH"
## MOB ------
load("./Data/MOB.dge.sceset.RData") 
spatial_count <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
spatial_count = t(spatial_count)
spatial_location = matrix(as.numeric(unlist(strsplit(colnames(spatial_count),'x'))),ncol=2,byrow=TRUE)
rownames(spatial_location) = colnames(spatial_count)
spatial_location = as.data.frame(spatial_location)
colnames(spatial_location) = c('x','y')

ct.varname = "cellType"
sample.varname = "sampleInfo"
ct.select = as.character(unique(pData(eset)$cellType))

CARD_obj = createCARDObject(
  sc_count = exprs(eset),
  sc_meta = pData(eset),
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'ST_MOB_ref.rds')

## HBC -------
spatial_count <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
spatial_count = t(spatial_count)
spatial_location = matrix(as.numeric(unlist(strsplit(colnames(spatial_count),'x'))),ncol=2,byrow=TRUE)
rownames(spatial_location) = colnames(spatial_count)
spatial_location = as.data.frame(spatial_location)
colnames(spatial_location) = c('x','y')

load(paste0(DATAPATH,"Azizi.BC.sceset"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)
print(unique(colData(eset)$sampleID))

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'ST_HBC_ref.rds')

## Slideseq ---------
load(paste0(DATAPATH,"/mouseCerebellum.sceset"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)

spatial_count <- readRDS(paste0(datapath,"Slide-seq.rds"))
spatial_count = spatial_count[which(rowSums(spatial_count)>0),which(colSums(spatial_count)>0)]

spatial_location = read.csv(paste0(datapath,"Slideseq_location.csv"))
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='barcodes')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
colnames(spatial_location) = c('x','y')

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'Slideseq_ref.rds')

## MC(C)
## 10x_HH -------- 
spatial_count<-Read10X_h5(paste0(datapath,"./HumanHeart10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_HH_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')

load(paste0(DATAPATH,"Human_Developmental_heart.sceset"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)
print(unique(colData(eset)$sampleID))

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_HH_ref.rds')

## 10x_MK------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseKidney10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_MK_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')

load(paste0(DATAPATH,"/mouseKidney.sceset"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_MK_ref.rds')

## MB(C) --------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrain10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_MB(C)_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')

iseed = 7
load(paste0(DATAPATH,"/Moubrain_part_l5_all_sub.rep",iseed,".RData"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_MB(C)_ref.rds')

## 10x_MB(S-A) ---------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrainAnterior10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_MB(S-A)_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')


iseed = 7
load(paste0(DATAPATH,"/Moubrain_part_l5_all_sub.rep",iseed,".RData"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_MB(S-A)_ref.rds')

## 10x_MB(S-P) ---------
spatial_count<-Read10X_h5(paste0(datapath,"./MouseBrainPosterior10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_MB(S-P)_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')

iseed = 7
load(paste0(DATAPATH,"/Moubrain_part_l5_all_sub.rep",iseed,".RData"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_MB(S-P)_ref.rds')


## 10x_HBC ------
spatial_count<-Read10X_h5(paste0(datapath,"./BC10x.h5"))
spatial_location=read.csv(paste0(datapath,"./10x_HBC_location.csv"),header = F)
spatial_location = spatial_location%>% remove_rownames %>% column_to_rownames(var='V1')
spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),]
spatial_location = spatial_location[,c('V5','V6')]
colnames(spatial_location) = c('x','y')

load(paste0(DATAPATH,"Azizi.BC.sceset"))
ct.varname = "cellType"
sample.varname = "sampleID"
ct.select = unique(colData(eset)$cellType)
print(unique(colData(eset)$sampleID))

CARD_obj = createCARDObject(
  sc_count = counts(eset),
  sc_meta = colData(eset)[,c('cellType','sampleID')],
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = ct.varname,
  ct.select = ct.select,
  sample.varname = sample.varname,
  minCountGene = 100,
  minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_obj)
saveRDS(CARD_obj,'10x_HBC_ref.rds')
