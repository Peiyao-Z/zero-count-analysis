library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(Matrix)
library(tidyverse)
library(tidyr)
library(dplyr)
library(gee)
library(geepack)
library(pscl)

datapath = "datapath"

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
  
  #number_of_cores= 4
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

  # newDat$mom_vsP = tmp$mom_vsP
  # newDat$mom_vsZINB = tmp$mom_vsZINB
  # newDat = data.frame(mom_vsP = tmp$mom_vsP,
  #                      mom_vsZINB = tmp$mom_vsZINB)
  
  return(AIC_sum)
}

## Niche ---
#dat = read.table(("NiCHE_GSM2788364_AB1655.txt"))
dat = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"NiCHE-seq_sum.rds")



# seqFISH+
dat<-read.csv(paste0(datapath,"ob_counts.csv"))
dat<-t(dat)
newDat = readRDS(paste0("AIC_offset/","seqFISH+","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/seqFISH+_sum.rds")
rm(dat,AIC_sum)


## # #STARmap ------------
dat<-read.csv(paste0(datapath,"STARmap.csv"),header = F)
dat<-t(dat)
newDat = readRDS(paste0("AIC_offset/","STARmap","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/STARmap_sum.rds")
rm(dat,AIC_sum)


## Slide-seq ----
dat<-readRDS(paste0(datapath,"Slide-seq.rds"))
newDat = readRDS(paste0("AIC_offset/","Slide-seq","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/Slide-seq_sum.rds")
rm(dat,AIC_sum)



## tomo-seq -------
dat<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
newDat = readRDS(paste0("AIC_offset/","Tomo-seq","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/Tomo-seq_sum.rds")
rm(dat,AIC_sum)




##LCM -----------
dat<-read.csv(paste0(datapath,"LCM.tsv"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
newDat = readRDS(paste0("AIC_offset/","LCM","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/LCM_sum.rds")
rm(dat,AIC_sum)


##liver single cell ------
dat<-read.csv(paste0(datapath,"GSE108561_NPC_umitab.txt"),sep = "\t")
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
# newDat = readRDS(paste0("AIC_offset/","PC-seq","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/PC-seq_sum.rds")
rm(dat,AIC_sum)


#ST -----------
dat<-read.csv(paste0(datapath,"L8CN18_D1_stdata_aligned_counts_IDs_ALS.txt"),sep = "\t",row.names = NULL)
dat<-dat[!duplicated(dat[,1]),]
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/ST_MALS_sum.rds")
rm(dat,AIC_sum)


dat <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
newDat = readRDS(paste0("AIC_offset/","ST_HBC","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/ST_HBC_sum.rds")
rm(dat,AIC_sum)


dat <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
dat <- t(dat)
newDat = readRDS(paste0("AIC_offset/","ST_MOB","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/ST_MOB_sum.rds")
rm(dat,AIC_sum)

#AdultMouseBrain10x
#cell.names<-read.table("MouseBrain10x/barcodes.tsv",stringsAsFactors = F)
dat<-readRDS(paste0(datapath,"MouseBrain10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_MB(C)","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_MB(C)_sum.rds")
rm(dat,AIC_sum)


#BC
dat<-readRDS(paste0(datapath,"BC10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_HBC","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_HBC_sum.rds")
rm(dat,AIC_sum)


#HumanHeart10x
dat<-readRDS(paste0(datapath,"HumanHeart10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_HH","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_HH_sum.rds")
rm(dat,AIC_sum)

#HumanLymph
dat<-readRDS(paste0(datapath,"HumanLymph10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_HL","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_HL_sum.rds")
rm(dat,AIC_sum)

#MouseBrainAnterior
dat<-readRDS(paste0(datapath,"MouseBrainAnterior10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_MB(S-A)","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_MB(S-A)_sum.rds")
rm(dat,AIC_sum)


#MouseBrainPosterior
dat<-readRDS(paste0(datapath,"MouseBrainPosterior10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_MB(S-P)","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_MB(S-P)_sum.rds")
rm(dat,AIC_sum)


#MouseKidney
dat<-readRDS(paste0(datapath,"MouseKidney10x_count.rds"))
newDat = readRDS(paste0("AIC_offset/","10x_MK","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/10x_MK_sum.rds")
rm(dat,AIC_sum)


# MERFISH
#setwd("D:/droplet/MERFISH")
dat<-read.csv(paste0(datapath,"MERFISH.csv"))
dat<-dat %>% remove_rownames %>% column_to_rownames(var=names(dat)[1])
dat<-t(dat)
newDat = readRDS(paste0("AIC_offset/","MERFISH","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/MERFISH_sum.rds")
rm(dat,AIC_sum)


##Slide-seqV2-----
dat<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
newDat = readRDS(paste0("AIC_offset/","Slide-seqV2","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/Slide-seqV2_sum.rds")
rm(dat,AIC_sum)

##seqFISH -----
library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==43)]
newDat = readRDS(paste0("AIC_offset/","seqFISH","_sum.rds"))
AIC_sum = total(dat)
saveRDS(AIC_sum,"AIC_offset/seqFISH_sum.rds")
rm(dat,AIC_sum)


#HDST-------
load(paste0(datapath,"CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds"))
dat<-sp_count
rm(sp_count)
newDat = readRDS(paste0("AIC_offset/","HDST","_sum.rds"))
st<-Sys.time()
AIC_sum = total(dat)
Sys.time()-st
saveRDS(AIC_sum,"AIC_offset/HDST_sum.rds")
rm(dat,AIC_sum)



