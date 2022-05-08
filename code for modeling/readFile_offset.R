library(parallel)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(MASS)
library(foreach)
library(doParallel)
library(Matrix)
library(Seurat)

library(pscl)
library(tidyverse)
library(vcdExtra)
library(lmtest)

datapath = "DEFINE_PATH_HERE"


### If we switch to source("function.R"), then it will be modeling count data without library size as offset/
source("function_offset.R")

##seqFISH -----
library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==34)]
source("gene-wise_offset.R")
saveRDS(sum_stat,"./seqFISH_sum.rds")


#STARmap-----
dat<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
dat<-t(dat)
source("gene-wise_offset.R")
saveRDS(sum_stat,"./STARmap_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#LCM-----
dat<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
source("gene-wise_offset.R")
saveRDS(sum_stat,"./LCM-seq_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#pcRNAseq-----
dat<-read.csv(paste0(datapath,"liver single cell.txt"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
source("gene-wise_offset.R")
saveRDS(sum_stat,"./pcRNAseq_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#ST
dat <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
source("gene-wise_offset.R")
saveRDS(sum_stat,"./ST_HBC_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

dat <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
source("gene-wise_offset.R")
saveRDS(sum_stat,"./ST_MOB_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#AdultMouseBrain10x
dat<-Read10X_h5(paste0(datapath,"./MouseBrain10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_MB(C)_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#BC
dat<-Read10X_h5(paste0(datapath,"./BC10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_HBC_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#HumanHeart10x
dat<-Read10X_h5(paste0(datapath,"./HumanHeart10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_HH_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#HumanLymph
dat<-Read10X_h5(paste0(datapath,"./HumanLymph10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_HL_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#MouseBrainAnterior
dat<-Read10X_h5(paste0(datapath,"./MouseBrainAnterior10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_MB(S-A)_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#MouseBrainPosterior
dat<-Read10X_h5(paste0(datapath,"./MouseBrainPosterior10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_MB(S-P)_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#MouseKidney
dat<-Read10X_h5(paste0(datapath,"./MouseKidney10x_count.h5"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./10x_MK_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#tomo-seq
dat<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
source("gene-wise_offset.R")
saveRDS(sum_stat,"./Tomo-seq_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#MERFISH
dat<-read.csv(paste0(datapath,"MERFISH.csv"))
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat<-t(dat)
source("gene-wise_offset.R")
saveRDS(sum_stat,"./MERFISH_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)


#seqFISH+
dat<-read.csv(paste0(datapath,"ob_counts.csv"))
dat<-t(dat)
source("gene-wise_offset.R")
saveRDS(sum_stat,"./seqFISH_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#NICHE-seq ------
dat = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"new/NiCHE-seq_sum.rds")


#Slide-seq ----
dat<-readRDS(paste0(datapath,"Slide-seq.rds"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./Slide-seq_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#Slide-seqV2-----
dat<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
source("gene-wise_offset.R")
saveRDS(sum_stat,"./Slide-seqV2_sum.rds")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

#HDST
load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
dat<-sp_count
rm(sp_count)
source("gene-wise_offset.R")
saveRDS(sum_stat,"new/HDST_sum.csv")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)