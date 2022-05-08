library(cowplot)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
summaryAICNB<-function(dat){
  sum_AIC<-data.frame(Var1=as.factor(c(1,2,3,4)))
  tol<-round(prop.table(table(dat[,c('final_chooseAIC')]))*100,2)
  tol<-as.data.frame(tol)
  sum_AIC<-left_join(sum_AIC,tol,by='Var1')
  inflate<-which(dat[,c('final_NBvsZINB')]<(0.05/nrow(dat)))
  inf<-round(prop.table(table(dat[inflate,c('final_chooseAIC')]))*100,2)
  inf<-as.data.frame(inf)
  if(nrow(inf)==0){
    inf<-data.frame(Var1 = as.factor(c(1,2,3,4)),Freq = factor(rep(NA,4)))
    levels(inf$Freq) <- c(levels(inf$Freq), 0)
    inf$Freq[is.na(inf$Freq)] <-0
    inf$Freq<-as.numeric(inf$Freq)-1
  }else{
  }
  tryCatch({
    sum_AIC<-left_join(sum_AIC,inf,by='Var1')
  }, error = function(error_condition) {
  })
  nonInflate<-which(dat[,c('final_NBvsZINB')]>=(0.05/nrow(dat)))
  nonInf<-round(prop.table(table(dat[nonInflate,c('final_chooseAIC')]))*100,2)
  nonInf<-as.data.frame(nonInf)
  if(nrow(nonInf)==0){
    nonInf<-data.frame(Var1 =as.factor(c(1,2,3,4)),Freq = factor(rep(NA,4)))
    levels(nonInf$Freq) <- c(levels(nonInf$Freq), 0)
    nonInf$Freq[is.na(nonInf$Freq)] <-0
  }else{
  }
  tryCatch({
    sum_AIC<-left_join(sum_AIC,nonInf,by='Var1')
  }, error = function(error_condition) {
  })
  colnames(sum_AIC)[1:2]<-c('model',"total")
  tryCatch({
    colnames(sum_AIC)[3]<-'inflated'
  }, error = function(error_condition) {
    
  })
  tryCatch({
    colnames(sum_AIC)[4]<-'non-inflated'
  }, error = function(error_condition) {
    
  })
  sum_AIC[is.na(sum_AIC)]<-0
  sum_AIC<-gather(sum_AIC,model,percentl)
  sum_AIC<-as.numeric(sum_AIC$percentl)
  #sum_AIC<-sum_AIC[,-which(colnames(sum_AIC)%in%'model')]
  
  if(length(sum_AIC)==12){
    return(sum_AIC)
  }else{
    #return(c(sum_AIC,rep(NA,(12-length(sum_AIC)))))
  }
  rm(tol,inf,nonInf,inflate,nonInflate)
}

list<-c("MERFISH","seqFISH","seqFISH+",
        "LCM-seq","NICHE-seq",'pcRNAseq',"STARmap","Tomo-seq",
        "HDST",'Slide-seq','Slide-seqV2',
        '10x_HBC','10x_HH','10x_HL','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')

## var vs mean & zero proportion vs mean (S1 & S2)-----

list<-c("MERFISH","seqFISH","STARmap",
        "LCM-seq","NICHE-seq","Tomo-seq",'pcRNAseq',
        "HDST",'Slide-seq',
        '10x_HBC','10x_HH','10x_HL','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')


zz<-list()
for (i in list){
  dat<-readRDS(paste0("./new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  dat$zeros_p<-exp(-dat$obs_mean)
  
  tmp = coef(nls(obs_zerofrac~(1+tmp*obs_mean)^(-1/tmp), data = dat,start = list(tmp = 1),trace = F,control = nls.control(warnOnly=T)))
  dat$zeros_nb<- (1+tmp*dat$obs_mean)^(-1/tmp)
  
  zz[[i]]<-ggplot(dat, mapping = aes(x = log(obs_mean), y = obs_zerofrac))+geom_point()+
    theme_classic() +labs(x='',y='',title = i)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line =element_line(colour = "black"),
          axis.title = element_blank(),
          axis.text = element_text(size = 12,face = "bold"),
          title = element_text(size=16))+
    geom_line(aes(x = log(obs_mean), y = zeros_nb), color = "red",size=1.3)+
    geom_line(aes(x = log(obs_mean), y = zeros_p), color = "blue",size=1.3)+
    theme(axis.title.y = element_text(size = 14,face = "bold"),axis.text = element_blank())
}

plt = ggarrange(zz[[1]],zz[[2]],zz[[3]],zz[[4]],zz[[5]],zz[[6]],zz[[7]],zz[[8]],zz[[9]],zz[[10]],zz[[11]],zz[[12]],zz[[13]],zz[[14]],zz[[15]],zz[[16]],zz[[17]])
plt = annotate_figure(plt,left = text_grob("Zero Proportion",face = "bold", size = 20,rot =90),
                      bottom = text_grob("Gene Mean",face = "bold", size = 20))

png('figS1.png',width = 30,height = 20,res = 2000,units = 'cm')
plot(plt)
dev.off()


pp<-list()
for (i in list){
  dat<-readRDS(paste0("./new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  tmp = coef(nls(obs_var ~ obs_mean + obs_mean^2/th, data = dat,start = list(th = 1),trace = F,control = nls.control(warnOnly=T)))
  dat$newVar<-dat$obs_mean+(1/tmp)*(dat$obs_mean)^2
  pp[[i]]<-ggplot(dat, mapping = aes(x = log(obs_mean), y = log(obs_var)))+geom_point()+
    theme_classic() +labs(x='',y='',title = i)+
    geom_abline(intercept=0,slope=1,color = "blue",size = 1.3)+
    geom_line(aes(x=log(obs_mean),y=log(newVar)),color = "red",size=1.3)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line =element_line(colour = "black"),
          axis.title = element_blank(),
          title = element_text(size=16),
          axis.title.y = element_text(size = 14,face = "bold"),axis.text = element_blank())
}

plt = ggarrange(pp[[1]],pp[[2]],pp[[3]],pp[[4]],pp[[5]],pp[[6]],pp[[7]],pp[[8]],pp[[9]],pp[[10]],pp[[11]],pp[[12]],pp[[13]],pp[[14]],pp[[15]],pp[[16]],pp[[17]])
plt = annotate_figure(plt,left = text_grob("Gene Variance",face = "bold", size = 20,rot =90),
                      bottom = text_grob("Gene Mean",face = "bold", size = 20))

png('figS2.png',width = 30,height = 20,res = 200,units = 'cm')
plot(plt)
dev.off()


## add cell type proportion as covariate S3 -----
list<-c('Slide-seq','10x_HBC','10x_HH','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')
pp=NULL
P_ratio = NULL  
NB_ratio = NULL  
ZIP_ratio = NULL  
ZINB_ratio = NULL 
num = NULL
for (i in list) {
  dat = readRDS(paste0("new2/",i,"_sum.rds"))
  summ = readRDS(paste0("CARD/",i,"_ref_sum.rds"))
  test = readRDS(paste0("CARD/",i,"_ref_test.rds"))
  ref = readRDS(paste0("CARD/",i,"_ref.rds"))
  
  pp = c(pp,dim(ref@Proportion_CARD)[1])
  num = c(num, dim(ref@Proportion_CARD)[2])
  summ<-summ[which(dat$obs_zerofrac<1),]
  test<-test[which(dat$obs_zerofrac<1),]
  dat<-dat[which(dat$obs_zerofrac<1),]
  
  P_ratio = c(P_ratio, length(which(summ$chooseAIC==1)) / length(which(dat$chooseAIC==1)))
  NB_ratio = c(NB_ratio, length(which(summ$chooseAIC==2)) / length(which(dat$chooseAIC==2)))
  ZIP_ratio = c(ZIP_ratio, length(which(summ$chooseAIC==3)) / length(which(dat$chooseAIC==3)))
  ZINB_ratio = c(ZINB_ratio, length(which(summ$chooseAIC==4)) / length(which(dat$chooseAIC==4)))
  
}

res1 = data.frame(Poisson = P_ratio, NB=NB_ratio,
                  ZIP=ZIP_ratio, ZINB=ZINB_ratio, tech=list)
res1 = gather(res1,model,val,-tech)
res1$model = factor(res1$model,levels = c('Poisson','NB','ZIP','ZINB'))
res1$tech = factor(res1$tech,levels=list)
fourA = ggplot(data = res1, mapping = aes(x = tech, y = (val), color = model,group = model))+ geom_point(size = 3,stat="identity")+
  geom_line(size=1.8)+#ylim(-3,3)+
  labs(y='Ratio of the preferred models \n with and without deconvolution',x='',title = "")+
  theme_bw()+#ylim(0,64)+
  geom_hline(aes(yintercept=c(1)),color = "bisque3",size = 1,key_glyph = "path",show.legend=FALSE)+
  theme(strip.text.x = element_text(size = 18),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.title=element_text(size=22),
        axis.text.y = element_text(size = 22,angle = 90,hjust=0.5),
        axis.text.x = element_text(size = 22,angle = 90,vjust=0.5),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'top')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(trans = 'log',
                     labels = scales::number_format(accuracy = 0.01L))



PvsNB = NULL
PvsZIP = NULL
NBvsZINB = NULL
ZIPvsZINB = NULL
for (i in list) {
  dat = readRDS(paste0("new2/",i,"_sum.rds"))
  summ = readRDS(paste0("CARD/",i,"_ref_sum.rds"))
  test = readRDS(paste0("CARD/",i,"_ref_test.rds"))
  ref = readRDS(paste0("CARD/",i,"_ref.rds"))
  
  pp = c(pp,dim(ref@Proportion_CARD)[1])
  num = c(num, dim(ref@Proportion_CARD)[2])
  summ<-summ[which(dat$obs_zerofrac<1),]
  test<-test[which(dat$obs_zerofrac<1),]
  dat<-dat[which(dat$obs_zerofrac<1),]
  
  PvsNB = c(PvsNB, length(which(test$PvsNB<(0.05/nrow(dat)))) / length(which(dat$PvsNB<(0.05/nrow(dat)))))
  PvsZIP = c(PvsZIP, length(which(test$PvsZIP<(0.05/nrow(dat)))) / length(which(dat$PvsZIP<(0.05/nrow(dat)))))
  NBvsZINB = c(NBvsZINB, length(which(test$NBvsZINB<(0.05/nrow(dat)))) / length(which(dat$NBvsZINB<(0.05/nrow(dat)))))
  ZIPvsZINB = c(ZIPvsZINB, length(which(test$ZIPvsZINB<(0.05/nrow(dat)))) / length(which(dat$ZIPvsZINB<(0.05/nrow(dat)))))
}

res = data.frame(PvsNB = PvsNB, PvsZIP=PvsZIP,
                 NBvsZINB=NBvsZINB, ZIPvsZINB=ZIPvsZINB, tech=list)
res = gather(res,model,val,-tech)
res$model = factor(res$model,levels = c('PvsZIP','PvsNB','NBvsZINB','ZIPvsZINB'))
res$tech = factor(res$tech,levels=list)

fourB = ggplot(data = res, mapping = aes(x = tech, y = (val), color = model,group = model))+ geom_point(size = 3,stat="identity")+
  geom_line(size=1.8)+#ylim(-3,3)+
  labs(y='Ratio of the significant genes \n with and without deconvolution',x='',title = "")+
  theme_bw()+
  geom_hline(aes(yintercept=c(1)),color = "bisque3",size = 1,key_glyph = "path",show.legend=FALSE)+
  theme(strip.text.x = element_text(size = 18),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.title=element_text(size=22),
        axis.text.y = element_text(size = 22,angle = 90,hjust=0.5),
        axis.text.x = element_text(size = 22,angle = 90,vjust=0.5),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'top')+
  # guides(colour = guide_legend(order = 1), 
  #        shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(trans = 'log',
                     labels = scales::number_format(accuracy = 0.01L))


plot4 = ggarrange(fourA,fourB, ncol = 1, nrow = 2)

ggsave("plot/S3.png",plot = plot4,width = 35,height = 32,
       units = 'cm', dpi = 1500)


## with and without AIC difference S4 -----
list<-c("seqFISH+",
        "LCM-seq","NICHE-seq",'pcRNAseq',"STARmap","Tomo-seq",
        "HDST",'Slide-seq','Slide-seqV2',
        '10x_HBC','10x_HH','10x_HL','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')

## no seqFISH, MERFISH
P_offset = NULL
P_no_offset = NULL
NB_offset = NULL
NB_no_offset = NULL
ZIP_offset = NULL
ZIP_no_offset = NULL
ZINB_offset = NULL
ZINB_no_offset = NULL

for (i in list) {
  dat1<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat1<-dat1[which(dat1$obs_zerofrac<1),]
  P_offset = c(P_offset,mean(dat1$P_aic,na.rm = T))
  NB_offset = c(NB_offset,mean(dat1$final_NB_aic,na.rm = T))
  ZIP_offset = c(ZIP_offset,mean(dat1$ZIP_aic,na.rm = T))
  ZINB_offset = c(ZINB_offset,mean(dat1$ZINB_aic,na.rm = T))
  dat2 = readRDS(paste0("no_offset/",i,"_sum.rds"))
  dat2<-dat2[which(dat2$obs_zerofrac<1),]
  P_no_offset = c(P_no_offset,mean(dat2$P_aic,na.rm = T))
  NB_no_offset = c(NB_no_offset,mean(dat2$NB_aic,na.rm = T))
  ZIP_no_offset = c(ZIP_no_offset,mean(dat2$ZIP_aic,na.rm = T))
  ZINB_no_offset = c(ZINB_no_offset,mean(dat2$ZINB_aic,na.rm = T))
}
tmp = data.frame(P_offset=P_offset,P_no_offset = P_no_offset,technology = list)
tmp$technology = factor(tmp$technology,levels = list)
aa1 = ggplot(tmp,aes(x=P_no_offset,y=P_offset,color=technology))+geom_point(size=5)+
  geom_abline(intercept=0,slope=1,color = "black",size = 1)+
  theme_bw()+ylim(0,6000)+xlim(0,6000)+
  theme(strip.text = element_blank(),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_blank(),
        axis.title=element_text(size=18,face = "bold"),
        axis.ticks = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22,face = "bold"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

tmp = data.frame(NB_offset=NB_offset,NB_no_offset = NB_no_offset,technology = list)
tmp$technology = factor(tmp$technology,levels = list)
aa2 = ggplot(tmp,aes(x=NB_no_offset,y=NB_offset,color=technology))+geom_point(size=5)+
  geom_abline(intercept=0,slope=1,color = "black",size = 1)+
  theme_bw()+ylim(0,6000)+xlim(0,6000)+
  theme(strip.text = element_blank(),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_blank(),
        axis.title=element_text(size=18,face = "bold"),
        axis.ticks = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22,face = "bold"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

tmp = data.frame(ZIP_offset=ZIP_offset,ZIP_no_offset = ZIP_no_offset,technology = list)
tmp$technology = factor(tmp$technology,levels = list)
aa3 = ggplot(tmp,aes(x=ZIP_no_offset,y=ZIP_offset,color=technology))+geom_point(size=5)+
  geom_abline(intercept=0,slope=1,color = "black",size = 1)+
  theme_bw()+ylim(0,6000)+xlim(0,6000)+
  theme(strip.text = element_blank(),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_blank(),
        axis.title=element_text(size=18,face = "bold"),
        axis.ticks = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22,face = "bold"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

tmp = data.frame(ZINB_offset=ZINB_offset,ZINB_no_offset = ZINB_no_offset,technology = list)
tmp$technology = factor(tmp$technology,levels = list)
aa4 = ggplot(tmp,aes(x=ZINB_no_offset,y=ZINB_offset,color=technology))+geom_point(size=5)+
  geom_abline(intercept=0,slope=1,color = "black",size = 1)+
  theme_bw()+ylim(0,6000)+xlim(0,6000)+
  theme(strip.text = element_blank(),legend.title= element_blank(),
        legend.text=element_text(size=20),
        axis.text = element_blank(),
        axis.title=element_text(size=18,face = "bold"),
        axis.ticks = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 22,face = "bold"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

plt4 = ggarrange(aa1,aa2,aa3,aa4,common.legend = T,legend  = "bottom")
png('plot/suppS4.png',width = 30,height = 20,units = 'cm',res=2000)
plot(plt4)
dev.off()

## celltype plots (S5)------
library(Seurat)
list<-c("MERFISH","seqFISH",
        "NICHE-seq",'pcRNAseq',"Tomo-seq",
        'Slide-seqV2',
        '10x_HBC','10x_HH','10x_HL','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')


pl = list()

# datapath = "DATAPATH"

clusterPC<-function(tmp,npcs = 15,resolution =0.5){
  colnames(tmp)<-c(1:ncol(tmp))
  tmp <- CreateSeuratObject(counts = tmp)
  tmp <- NormalizeData(tmp)
  all.genes <- rownames(tmp)
  tmp <- ScaleData(tmp, features = all.genes)
  tmp <- FindVariableFeatures(object = tmp)
  tmp <- RunPCA(tmp, features = VariableFeatures(object = tmp),npcs = npcs)
  tmp<-FindNeighbors(tmp,dims=1:10)
  tmp<-FindClusters(tmp,resolution = resolution )
  tmp <- RunUMAP(tmp, dims = 1:10)
  return(tmp)
}

## seqFISH -----
library(openxlsx)
dat = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=1,colNames=F,rowNames = T)
tmp = read.xlsx(paste0(datapath,"mmc6.xlsx"),sheet=2,colNames=F,rowNames = T)
dat = dat[,which(tmp[1,]==43)]
data = dat
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["seqFISH"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "seqFISH")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


## MERFISH -----
data<-read.csv(paste0(datapath,"MERFISH.csv"))
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
data<-t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["MERFISH"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "MERFISH")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))

## ST----
data <- read.table(paste0(datapath,"Layer2_BC_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["ST_HBC"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "ST_HBC")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


data <- read.table(paste0(datapath,"Rep11_MOB_count_matrix-1.tsv"),sep = "")
data <- t(data)
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["ST_MOB"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "ST_MOB")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


#NICHE-seq ------
dat = read.table(paste0(datapath,"NiCHE_GSM2788364_AB1655.txt.gz"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["NICHE-seq"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "ST_HBC")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))
## Tomo-seq ----
data<-read.table(paste0(datapath,"GSM3148567_Male.animal.1.ReadCounts.tsv"),stringsAsFactors = F,header = T)
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["Tomo-seq"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "ST_HBC")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))
## 10x-----------
#BC
data<-Read10X_h5(paste0(datapath,"10x/BC10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_HBC"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_HBC")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


data<-Read10X_h5(paste0(datapath,"10x/HumanHeart10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_HH"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_HH")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


data<-Read10X_h5(paste0(datapath,"10x/HumanLymph10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_HL"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_HL")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))



data<-Read10X_h5(paste0(datapath,"10x/MouseKidney10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_MK"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_MK")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))

data<-Read10X_h5(paste0(datapath,"10x/MouseBrain10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_MB(C)"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_MB(C)")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))

data<-Read10X_h5(paste0(datapath,"10x/MouseBrainAnterior10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_MB(S-A)"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_MB(S-A)")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


data<-Read10X_h5(paste0(datapath,"10x/MouseBrainPosterior10x.h5"))
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["10x_MB(S-P)"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "10x_MB(S-P)")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


## Paired_cell_sequencing ----
data<-read.csv(paste0(datapath,"liver single cell.txt"),sep = "\t")
data<-data %>% remove_rownames %>% column_to_rownames(var=names(data)[1])
tmp<-clusterPC(tmp = data,npcs = 15,resolution = 0.5)
pl[["PC-seq"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "pcRNAseq")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))


## Slide-seq2 ----- 
data<-readRDS(paste0(datapath,"Slide-seqV2.rds"))
tmp<-clusterPC(tmp = data)
pl[["Slide-seqV2"]] = DimPlot(object = tmp, reduction = "umap")+labs(title = "Slide-seqV2")+
  theme(panel.grid =element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        title = element_text(size=18,face = "bold"))

UMAP = plot_grid(pl[[1]],pl[[2]],pl[[3]],pl[[4]],
          pl[[5]],pl[[6]],pl[[7]],pl[[8]],
          pl[[9]],pl[[10]],pl[[11]],pl[[12]],
          pl[[13]],pl[[14]],pl[[15]])
png('plot/suppS5.png',width = 30,height = 20,units = 'cm',res=2000)
plot(UMAP)
dev.off()
