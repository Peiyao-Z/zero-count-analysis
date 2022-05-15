
library(tidyverse)
library(Seurat)

zero = NULL
readDepth = NULL

datapath = "/net/mulan/disk2/peiyao/"


library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==34)]
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
dat<-t(dat)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))


dat<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-read.csv(paste0(datapath,"liver single cell.txt"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/MouseBrain10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/BC10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/HumanHeart10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/HumanLymph10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/MouseBrainAnterior10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/MouseBrainPosterior10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-Read10X_h5(paste0(datapath,"10x/MouseKidney10x.h5"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-read.csv(paste0(datapath,"MERFISH.csv"))
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat<-t(dat)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-read.csv(paste0(datapath,"ob_counts.csv"))
dat<-t(dat)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, (sum(dat) / ncol(dat)))

dat<-readRDS(paste0(datapath,"Slide-seq.rds"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

dat<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,sum(dat==0)/(nrow(dat)*ncol(dat)))
readDepth = c(readDepth, sum(dat) / ncol(dat))

load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
dat<-sp_count
rm(sp_count)
dat = dat[which(rowSums(dat)>0),which(colSums(dat)>0)]
zero = c(zero,0.999599)
# sum(dat==0) 3616821526
# dim(dat): 19950*181367
readDepth = c(readDepth, sum(dat) / ncol(dat))


list = c('seqFISH','STARmap','LCM-seq','pcRNAseq','ST_HBC','ST_MOB',
         '10x_MB(C)','10x_HBC','10x_HH','10x_HL','10x_MB(S-A)','10x_MB(S-P)','10x_MK',
         'Tomo-seq','MERFISH','seqFISH+','NICHE-seq','Slide-seq','Slide-seqV2','HDST')

tmp = data.frame(model = list, zeroProp = zero, readDepth=readDepth)
saveRDS(tmp,'readDepth_and_zeroProp.rds')

cor.test(zero[,2], zero[,3],method = "spearman")



