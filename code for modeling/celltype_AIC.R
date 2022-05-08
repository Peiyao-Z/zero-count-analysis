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

datapath = "DATA_PATH"

source("function_offset.R")
total = function(dat){
  
  tt = colSums(dat)
  
  warn = function(arr,tt){
    tmp = NULL
    arr = as.numeric(arr)
    tryCatch({
      nb = glm.nb(arr ~ 1+offset(log(tt)))
      tmp = 0
    },warning = function(w) {
      tmp <<- 1
    },error = function(e) {
      tmp <<- 1
    })
    return(tmp)
  }
  #
  err = function(arr,tt){
    tmp = NULL
    arr = as.numeric(arr)
    tryCatch({
      nb = glm.nb(arr ~ 1+offset(log(tt)))
      tmp = 0
    },error = function(e) {
      tmp <<- 1
    })
    return(tmp)
  }
  
  
  lik = function(arr,tt){
    tryCatch({
      arr=as.numeric(arr)
      #bb = geeglm(arr~1,family = poisson(link = "log"),id=c(1:length(arr)),offset = log(tt))
      tmp <- gee(arr~1+offset(log(tt)), id=c(1:length(arr)), family = poisson(link = "log"))
      #mu = bb$fitted.values
      mu = exp(coef(tmp)+log(tt))
      yy = mu*tmp$scale
      model<- lm(yy ~  1* mu + I(mu^2) + 0 )
      th=1/coef(model)
      #th = coef(nls(yy ~ mu + mu^2/th, start = list(th = 1),trace = F,control = nls.control(warnOnly=T)))
      loglik <- function(th, mu, y) sum( (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                            log(mu + (y == 0)) - (th + y) * log(th + mu)))
      AIC = 2*2-2*loglik(th = th, mu = mu, y=arr)
      
      nb = loglik(th = th, mu = mu, y=arr)
      po = logLik(glm(arr ~ offset(log(tt)),family = "poisson"))
      znb = logLik(zeroinfl(arr ~ 1+offset(log(tt))|1,dist = "negbin"))
      
      return(data.frame(mu = mean(mu),AIC=AIC,dispersion = th,
                        mom_vsP = pchisq(-2*(po-nb),df=1),mom_vsZINB = pchisq(-2*(nb-znb),df=1)))
    },error=function(e){
      return(data.frame(mu = NA,AIC = NA,dispersion = NA,
                        mom_vsP = NA, mom_vsZINB = NA))
    })
  }
  
  compare = function(arr,tt){
    tryCatch({
      arr=as.numeric(arr)
      tmp <- gee(arr~1+offset(log(tt)), id=c(1:length(arr)), family = poisson(link = "log"))
      mu = exp(coef(tmp)+log(tt))
      yy = mu*tmp$scale
      model<- lm(yy ~  1* mu + I(mu^2) + 0 )
      th=1/coef(model)
      loglik <- function(th, mu, y) sum( (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                            log(mu + (y == 0)) - (th + y) * log(th + mu)))
      nb = loglik(th = th, mu = mu, y=arr)
      po = logLik(glm(arr ~ offset(log(tt)),family = "poisson"))
      znb = logLik(zeroinfl(arr ~ 1+offset(log(tt))|1,dist = "negbin"))
      
      pchisq(-2*(po-nb),df=1)
      pchisq(-2*(nb-znb),df=1)
      return(data.frame(mom_vsP = pchisq(-2*(po-nb),df=1),mom_vsZINB = pchisq(-2*(nb-znb),df=1)))
    },error=function(e){
      return(data.frame(mom_vsP = NA, mom_vsZINB = NA))
    })
  }
  
  number_of_cores <- parallel::detectCores() - 1
  clusters <- parallel::makeCluster(number_of_cores)
  doParallel::registerDoParallel(clusters)
  warn2 = foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS'))%dopar% warn(dat[i,],tt)
  error2 = foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS'))  %dopar%  err(dat[i,],tt)
  lik = foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS','geepack','gee','pscl'))  %dopar%  lik(dat[i,],tt)
  tmp = foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS','geepack','gee','pscl')) %dopar% compare(dat[i,],tt = tt)
  parallel::stopCluster(clusters)
  
  tt = as.numeric(warn2)
  tt = ifelse(as.numeric(error2)==1,NA,tt)
  
  AIC_sum = data.frame(warning = tt,
                       mu = lik$mu,
                       dispersion=lik$dispersion,
                       AIC=lik$AIC,
                       mom_vsP = tmp$mom_vsP,
                       mom_vsZINB = tmp$mom_vsZINB)
  
  return(AIC_sum)
}

clusterPC<-function(data){
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

## HDST ----
load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
data<-sp_count
rm(location,sp_count)
ann = read.csv(paste0(datapath,"annotation.csv"))
ann$new = paste0(ann$spot_x,"x",ann$spot_y)
group = unique(ann$poly.ID)
cellType_sum <- list()
for (cellType in group[3:length(group)]) {
  tt = ann$new[which(ann$poly.ID==cellType)]
  #tt<-cluster$bc[which(cluster$cell_type==cellType)]
  dat <- data[,colnames(data)%in%tt]
  cellType_sum[[cellType]] <- total(dat)
  print(group)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info,tt)
}
saveRDS(cellType_sum,"celltype_mom/HDST_celltype.rds")
rm(data,tmp,cellType_sum)


## Slide-seq2 ----
data<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
tmp<-clusterPC(data = data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))

cellType_sum <- list()

for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/Slide-seqV2_celltype.rds")
rm(data,tmp,cellType_sum)



# STARmap----
data<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
data<-t(data)
cluster<-read.csv(paste0(datapath,"class_labels.csv"))
clu<-data.frame(CellID = 0:(ncol(data)-1))
cluster<-left_join(clu,cluster,by="CellID")
colnames(data)<-cluster$ClusterName
group<-names(sort(table(cluster$ClusterName),decreasing = T))

cellType_sum = list()

for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/STARmap_celltype.rds")
rm(data,dat,tmp,cellType_sum)


# ##seqFISH -----
library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==34)]
data = dat
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
  print(group)
}
saveRDS(cellType_sum,"celltype_mom/seqFISH_celltype.rds")
rm(data,dat,tmp,cellType_sum)


## seqFISH+ -----
data<-read.csv(paste0(datapath,"ob_counts.csv"))
data<-t(data)
cluster<-read.csv(paste0(datapath,"OB_cell_type_annotations.csv"))
cellType_sum <- list()
for (cellType in 1:25) {
  #print(cellType)
  dat<-data[,which(cluster$louvain==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/seqFISH+_celltype.rds")
rm(data,dat,tmp,cellType_sum)

## ST ----
data <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(data = data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))

cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/ST_HBC_celltype.rds")
rm(data,dat,tmp,cellType_sum)

data <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(data = data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))

cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/ST_MOB_celltype.rds")
rm(data,dat,tmp,cellType_sum)



# ## tomo -----
data<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(data = data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))

cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/Tomo-seq_celltype.rds")
rm(data,dat,tmp,cellType_sum)



# Paired_cell_sequencing ----
data<-read.csv(paste0(datapath,"liver single cell.txt"),sep = "\t")
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(data = data)
cellType_sum <- list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <-total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/pcRNAseq_celltype.rds")
rm(data,dat,tmp,cellType_sum)


## LCM -----
data<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
cluster<-read.csv(paste0(datapath,"cell_zone_table.txt"),sep = "")
colnames(data)<-cluster$zone
#group<-levels(as.factor(cluster$zone))
group<-names(sort(table(cluster$zone),decreasing = T))
cellType_sum <- list()
#cellType_sum <- list()
for (cellType in group) {
  dat <- data[,colnames(data)%in%cellType]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
  print(group)
}
saveRDS(cellType_sum,"celltype_mom/LCM_celltype.rds")
rm(data,dat,tmp,cellType_sum)

## MERFISH -----
data<-read.csv(paste0(datapath,"MERFISH.csv"))
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
data<-t(data)
tmp<-clusterPC(data = data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/MERFISH_celltype.rds")
rm(data,dat,tmp,cellType_sum)


# ## Slide-seq -----
data<-readRDS(paste0(datapath,"Slide-seq.rds"))
cluster<-read.csv(paste0(datapath,"AnalogizerClusterAssignments.csv"),header = T)
#cluster<-cluster[,c("bc","cluster")]
group<-unique(cluster$x)
cellType_sum <- list()
for (cellType in 1:11) {
  tt<-cluster$Var1[cluster$x==cellType]
  dat<-data[,colnames(data)%in%tt]
  cellType_sum[[cellType]] <- total(dat)
  print(group)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info,tt)
}
saveRDS(cellType_sum,"celltype_mom/Slide-seq_celltype.rds")
rm(data,dat,tmp,cellType_sum)

## 10x-----------
#AdultMouseBrain10x
data<-Read10X_h5(paste0(datapath,"MouseBrain10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_MB(C)_celltype.rds")
rm(data,dat,tmp,cellType_sum)

#BC
data<-Read10X_h5(paste0(datapath,"BC10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_HBC_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-Read10X_h5(paste0(datapath,"HumanHeart10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_HH_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-Read10X_h5(paste0(datapath,"HumanLymph10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_HL_celltype.rds")
rm(data,dat,tmp,cellType_sum)



data<-Read10X_h5(paste0(datapath,"MouseKidney10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_MK_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-Read10X_h5(paste0(datapath,"MouseBrainAnterior10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_MB(S-A)_celltype.rds")
rm(data,dat,tmp,cellType_sum)


data<-Read10X_h5(paste0(datapath,"MouseBrainPosterior10x.h5"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/10x_MB(S-P)_celltype.rds")
rm(data,dat,tmp,cellType_sum)


##NICHE ---------
data = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
tmp<-clusterPC(data=data)
#colnames(data)<-Idents(tmp)
group<-levels(Idents(tmp))
cellType_sum = list()
for (cellType in group) {
  dat <- data[,which(Idents(tmp)==cellType)]
  cellType_sum[[cellType]] <- total(dat)
  rm(dat,sum_stat,po_info,nb_info,zpo_info,znb_info)
}
saveRDS(cellType_sum,"celltype_mom/NICHE-seq_celltype.rds")
rm(data,dat,tmp,cellType_sum)


