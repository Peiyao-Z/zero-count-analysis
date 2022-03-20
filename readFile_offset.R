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
# library(pscl,lib.loc = "/home/pyzhao/R/x86_64-pc-linux-gnu-library/3.6")
# library(tidyverse,lib.loc = "/home/pyzhao/R/x86_64-pc-linux-gnu-library/3.6")
# library(vcdExtra,lib.loc = "/home/pyzhao/R/x86_64-pc-linux-gnu-library/3.6")
# library(lmtest,lib.loc = "/home/pyzhao/R/x86_64-pc-linux-gnu-library/3.6")

library(pscl)
library(tidyverse)
library(vcdExtra)
library(lmtest)

datapath = "/net/mulan/disk2/peiyao/"

source("function_offset.R")

# ##NICHE ---------
# dat = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/NiCHE-seq_sum.rds")
# 
# ## 
# # load(paste0(datapath,'NicheDataLCM.rda'))
# # dat = NicheDataLCM
# # rm(NicheDataLCM)
# # source("gene-wise_offset.R")
# # saveRDS(sum_stat,"new/LCM_Niche_sum.rds")
# 
# # ##seqFISH -----
# library(openxlsx)
# dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
# tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
# dat = dat[,which(tmp[1,]==34)]
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/seqFISH_sum.rds")
# 
# 
# #STARmap
# #setwd("D:/droplet/STARmap/")
# dat<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
# dat<-t(dat)
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/STARmap_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# 
# # #LCM
# # #setwd("D:/droplet/LCM/")
# dat<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
# dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/LCM-seq_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# #liver single cell
# #setwd("D:/droplet/liver single cell/")
# dat<-read.csv(paste0(datapath,"GSE108561_NPC_umitab.txt"),sep = "\t")
# dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
# #source("D:/droplet/gene-wise_offset.R")
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/pcRNAseq_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# 
# # #ST
# # #setwd("D:/droplet/MOB - ST/")
# # dat<-dat<-read.csv("L8CN18_D1_stdata_aligned_counts_IDs_ALS.txt",sep = "\t",row.names = NULL)
# # dat<-dat[!duplicated(dat[,1]),]
# # dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
# # #source("D:/droplet/gene-wise_offset.R")
# # source("gene-wise_offset.R")
# # saveRDS(sum_stat,"new/ST_ALS_sum.csv")
# # rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# dat <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
# dat <- t(dat)
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/ST_HBC_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# dat <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
# dat <- t(dat)
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/ST_MOB_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# #AdultMouseBrain10x
# #cell.names<-read.table("MouseBrain10x/barcodes.tsv",stringsAsFactors = F)
# dat<-readRDS(paste0(datapath,"10x/MouseBrain10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_MB(C)_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# #BC
# dat<-readRDS(paste0(datapath,"10x/BC10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_HBC_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# 
# #HumanHeart10x
# dat<-readRDS(paste0(datapath,"10x/HumanHeart10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_HH_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# #HumanLymph
# dat<-readRDS(paste0(datapath,"10x/HumanLymph10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_HL_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# #MouseBrainAnterior
# dat<-readRDS(paste0(datapath,"10x/MouseBrainAnterior10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_MB(S-A)_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# 
# #MouseBrainPosterior
# dat<-readRDS(paste0(datapath,"10x/MouseBrainPosterior10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_MB(S-P)_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
# 
# 
# #MouseKidney
# dat<-readRDS(paste0(datapath,"10x/MouseKidney10x_count.rds"))
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/10x_MK_sum.rds")
# rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

## # tomo-seq
#dat<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
#dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
#source("gene-wise_offset.R")
#saveRDS(sum_stat,"new/Tomo-seq_sum.rds")
#rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
#
#
## # MERFISH
## #setwd("D:/droplet/MERFISH")
#dat<-read.csv(paste0(datapath,"MERFISH.csv"))
#dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
#dat<-t(dat)
#source("gene-wise_offset.R")
#saveRDS(sum_stat,"new/MERFISH_sum.rds")
#rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
#
#
## # seqFISH+
## #setwd("D:/droplet/seqFISH")
#dat<-read.csv(paste0(datapath,"ob_counts.csv"))
#dat<-t(dat)
#source("gene-wise_offset.R")
#saveRDS(sum_stat,"new/seqFISH_sum.rds")
#rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
#
### Slide-seq ----
##setwd("D:/droplet/Slideseq")
#dat<-readRDS(paste0(datapath,"Slide-seq.rds"))
#source("gene-wise_offset.R")
#saveRDS(sum_stat,"new/Slide-seq_sum.rds")
#rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

###Slide-seqV2-----
#dat<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
#source("gene-wise_offset.R")
#saveRDS(sum_stat,"new/Slide-seqV2_sum.rds")
#rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)

# dat = Read10X(data.dir = "/net/mulan/disk2/peiyao/seqScope/liverData/")
# source("gene-wise_offset.R")
# saveRDS(sum_stat,"new/seqScope.rds")

#HDST
load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
dat<-sp_count
rm(sp_count)
source("gene-wise_offset.R")
saveRDS(sum_stat,"new/HDST_sum.csv")
rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)