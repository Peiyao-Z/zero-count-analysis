rm(list=ls())

list<-c("MERFISH","seqFISH","seqFISH+","STARmap",
        "LCM-seq","NICHE-seq","Tomo-seq","pcRNAseq",
        "HDST",'Slide-seq','Slide-seqV2',
        '10x_HBC','10x_HH','10x_HL','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')


chooseAIC<-function(arr){
  if(sum(!is.na(arr))==0){
    return(NA)
  }
  else{
    t = as.numeric(which.min(arr))
    return(t)
  }
}


for (i in list) {
  dat1 = readRDS(paste0("new/",i,"_sum.rds"))
  dat2 = readRDS(paste0("AIC_offset/",i,"_sum.rds"))
  dat1$NB_warning = dat2$warning
  dat1$NB_mom_mean = dat2$mu
  dat1$NB_mom_aic = dat2$AIC
  dat1$NB_mom_dispersion = dat2$dispersion
  dat1$NB_mom_var = dat1$NB_mom_mean*(1+dat1$NB_mom_mean/dat1$NB_mom_dispersion)
  dat1$NB_mom_exp_zerofrac = (1+dat1$NB_mom_mean/dat1$NB_mom_dispersion)^(-dat1$NB_mom_dispersion)
  dat1$NB_mom_PvsNB = dat2$mom_vsP
  dat1$NB_mom_NBvsZINB = dat2$mom_vsZINB

  dat1$final_NB_aic = apply(dat1[,c('NB_warning','NB_aic','NB_mom_aic')],1,
                            function(arr) ifelse(arr[1]==1,min(arr[2],arr[3],na.rm = T),arr[2]))
  
  dat1$final_NB_mean = apply(dat1[,c('final_NB_aic','NB_aic','NB_mom_aic','NB_mean','NB_mom_mean')],1,
                             function(arr) ifelse(arr[1]==arr[2],arr[4],arr[5]))
  
  dat1$final_NB_var = apply(dat1[,c('final_NB_aic','NB_aic','NB_mom_aic','NB_var','NB_mom_var')],1,
                             function(arr) ifelse(arr[1]==arr[2],arr[4],arr[5]))
  
  dat1$final_NB_exp_zerofrac = apply(dat1[,c('final_NB_aic','NB_aic','NB_mom_aic','NB_exp_zerofrac','NB_mom_exp_zerofrac')],1,
                             function(arr) ifelse(arr[1]==arr[2],arr[4],arr[5]))

  dat1$final_chooseAIC = apply(dat1[,paste0(c('P','final_NB','ZIP','ZINB'),'_aic')],1,chooseAIC)
  dat1$final_PvsNB = apply(dat1[,c('final_NB_aic','NB_aic','NB_mom_aic','PvsNB','NB_mom_PvsNB')],1,
                           function(arr) ifelse(arr[1]==arr[2],arr[4],arr[5]))
  dat1$final_NBvsZINB = apply(dat1[,c('final_NB_aic','NB_aic','NB_mom_aic','NBvsZINB','NB_mom_NBvsZINB')],1,
                              function(arr) ifelse(arr[1]==arr[2],arr[4],arr[5]))
  saveRDS(dat1,paste0("new2/",i,"_sum.rds"))
  print(i)
}

