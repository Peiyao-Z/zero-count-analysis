library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(EnvStats)
library(ggpubr)
library(grid)
library(ggpointdensity)
library(MASS)
library(gtools)
library(uwot)
require(scales)
library(tidyverse)
library(RColorBrewer)

list<-c("MERFISH","seqFISH","seqFISH+",
        "LCM-seq","NICHE-seq",'pcRNAseq',"STARmap","Tomo-seq",
        "HDST",'Slide-seq','Slide-seqV2',
        '10x_HBC','10x_HH','10x_HL','10x_MB(C)','10x_MB(S-A)','10x_MB(S-P)','10x_MK','ST_HBC','ST_MOB')

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

tt = NULL
for(i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  model<- lm(obs_var ~  1* obs_mean + I(obs_mean^2) + 0, data =dat)
  phi=coef(model)
  tt = c(tt,phi)
}
#total$theta = tt
total = data.frame(model=list,theta = tt)
saveRDS(total,"./overdispersionParam.rds")


fourList<-c('seqFISH+','Slide-seqV2','10x_MB(C)')

## figure 1 ----
zero<-readRDS("readDepth_and_zeroProp.rds")
zero$model<-factor(zero$model,levels = list)
oneA = ggplot(data = zero, mapping = aes( x=model,y = zeroProp, fill = model))+ geom_col(position =  "stack")+
  labs(y='Zero Proportion',x="")+theme_bw()+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "darkgray",size = 1,key_glyph = "path",show.legend = F)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell"),color = "darkgray",size = 1,key_glyph = "path",show.legend = F)+
  theme(axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 90,vjust = 0,size=20,face = "bold",hjust = 0.5),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size=20,face = "bold"),
        # axis.text.y = element_text(size=20,face = "bold"),
        axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(legend.position = "top",
      legend.title=element_blank(),
      legend.text = element_text(size=20),
      plot.margin = unit(c(0,30,0,0), "pt"),
      legend.box.margin=margin(20,20,20,20))+
  guides(color=guide_legend(nrow=3))

ggsave('oneA.png',dpi = 2000,width = 7.53 , height = 5.72)
rm(tmp,zero)

tmp<-readRDS("./readDepth_and_zeroProp.rds")
tmp$model<-factor(tmp$model,levels = list)
# tmp = tmp[tmp$model!="Tomo-seq",]
tt=tmp[which(!(tmp$model%in%c("MERFISH","seqFISH",'Tomo-seq'))),]
oneB = ggplot()+ 
  geom_point(data = tmp, mapping = aes(x=readDepth,y = logit(zeroProp),colour=model),size = 7)+
  #geom_text(aes(label=model),hjust=0, vjust=1,size = 6)+
  labs(y='Zero Proportion',x="Read depth per location")+theme_bw()+
  xlim(0,42000)+
  theme(legend.text=element_text(size=20,color = 'black'),
        axis.title=element_text(size=20,color = 'black',face = 'bold'),
        axis.text = element_text(size=20,color = 'black'),
        title = element_text(size=20,color = 'black'),
        #legend.title = element_blank(),
        legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(labels=c("-4" = "0.02", "0" = "0.50","4" = "0.98", "8" = "1.00"))+
  geom_smooth(data = tt,mapping = aes(x=readDepth,y = logit(zeroProp)),color = "black",alpha=0.25)+
  theme(legend.position = "top",
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        plot.margin = unit(c(0,30,0,0), "pt"),
        legend.box.margin=margin(20,20,20,20))+
  guides(color=guide_legend(nrow=3))
rm(tmp)
# 
# ggplot(data = tmp,mapping = aes(x=readDepth,y = logit(zeroProp),colour=model))+ 
#   geom_point( size = 7)+
#   geom_text(data = tmp,aes(label=model),hjust=0, vjust=1,size = 6)+
#   labs(y='Zero Proportion',x="Read depth per location")+theme_bw()

ggsave('oneB.png',dpi = 2000,width = 7.53 , height = 5.72)

plt = ggarrange(oneA,oneB,common.legend = TRUE, legend="top")
ggsave('plot/oneAB.png',dpi = 1500,width = 15 , height = 5.72)



fourList<-c('seqFISH+','Slide-seqV2','10x_MB(C)')
meanZero<-list()
for(i in fourList[1:3]){
  dat<-readRDS(paste0("./new2/",i,"_sum.rds"))
  dat<-dat[complete.cases(dat),]
  dat<-dat[which(dat$obs_zerofrac<1),]
  dat$zeros_p<-exp(-dat$obs_mean)
  phi = coef(nls(obs_var ~ obs_mean + obs_mean^2/th, data = dat,start = list(th = 1),trace = F,control = nls.control(warnOnly=T)))
  dat$zeros_nb<- (1+phi*dat$obs_mean)^(-1/phi)
  
  meanZero[[i]]<-ggplot(dat, mapping = aes(x = log(obs_mean), y = obs_zerofrac))+geom_point()+
    theme_bw() +labs(x='Gene Mean',y='Zero Proportion',title = i)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line =element_line(colour = "black"),
          axis.title = element_blank(),
          axis.text = element_text(size = 16,face = "bold"),
          title = element_text(size=18))+
    geom_line(aes(x = log(obs_mean), y = zeros_nb), color = "red",size=1.3)+
    geom_line(aes(x = log(obs_mean), y = zeros_p), color = "blue",size=1.3)
  rm(tmp,dat)
}

oneC<-plot_grid(meanZero[[1]], meanZero[[2]], meanZero[[3]],ncol=3,align = 'v')
oneC = annotate_figure(oneC,left = text_grob("Zero Proportion",face = "bold", size = 20,rot =90),
                       bottom = text_grob("Gene Mean",face = "bold", size = 20))
#ggsave('oneC.png',dpi = 2000,width = 7.53 , height = 5.72)
png('plot/oneC.png',width = 7.53 , height = 5.72,res=2000)
plot(oneC)
dev.off()

mseVarP<-NULL
mseZeroP<-NULL
mseVarNB<-NULL
mseZeroNB<-NULL

chart = data.frame(row.names = c("id","mseVarP","mseVarNB","mseZeroP","mseZeroNB"))
total = readRDS("./overdispersionParam.rds")

for(i in list){
  dat<-readRDS(paste0("./new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  #dat$label<-ntile(dat$obs_mean,4)
  tmp = total[which(total$model==i),]
  dat$zeros_p<-exp(-dat$obs_mean)
  dat$zeros_nb<- (1+dat$obs_mean*tmp$theta)^(-1/tmp$theta)
  mseZeroP<-c(mseZeroP,mean((dat$zeros_p - dat$obs_zerofrac))^2)
  mseZeroNB<-c(mseZeroNB,mean((dat$zeros_nb - dat$obs_zerofrac))^2)
}

mseZero<-data.frame(poisson=mseZeroP,NB=mseZeroNB)
mseZero<-gather(mseZero,model,Val)
mseZero$model<-factor(mseZero$model,levels = c("poisson","NB"))
mseZero$tech<-factor(rep(list,2),levels = list)
oneD = ggplot(data = mseZero, mapping = aes(x = tech, y = (Val), colour = model,group = model))+ geom_line(size = 1.5)+geom_point(size=2.5)+
  labs(y='MSE of Zero Proportion',x='')+theme_bw()+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  theme(axis.title=element_text(size=20,face = "bold"),
        legend.position = c(0.85, 0.25),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20,face = "bold",angle = 90,vjust = 0.4),axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=20,face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=17, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  scale_colour_manual(values = c( "blue","red"))+
  scale_y_continuous(trans = 'log',
                     labels = (scientific))


ggsave('./oneD.png',dpi = 2000,width = 6.96 , height = 5.72)

plot_grid(oneA,oneC,oneB,oneD,labels = c('A', 'B','C','D'),nrow=2,align = 'v')


## plot2----
pp = NULL
for(i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  pp = c(pp,length(which(dat$obs_var>dat$obs_mean))/nrow(dat))
}
tmp = data.frame(proportion=pp)
tmp$model = factor(list,levels = list)
twoA = ggplot(data = tmp, mapping = aes(x = model, y = proportion,color=model))+ geom_point(size = 5)+#geom_line(size=1.3,group = 1)+
  labs(y='Genes with variance > mean (%)',x='')+theme_bw()+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  geom_vline(aes(xintercept=c(8.5),linetype = "non-UMI data"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  theme(axis.title=element_text(size=20,face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=20),
        plot.title = element_text (hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  theme(legend.position = "top",
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        plot.margin = unit(c(0,30,0,0), "pt"),
        legend.box.margin=margin(20,20,20,20))+
  guides(color=guide_legend(nrow=3))

ggsave('twoA.png',dpi = 2000,width = 7.3 , height = 5.72)


mmedian = list()
for(i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  mmedian[[i]] = dat$obs_var/dat$obs_mean
}

tt <- data.frame(x = unlist(mmedian), 
                 model = rep(list[1:length(mmedian)],times = sapply(mmedian,length)))
tt$model = factor(tt$model,levels = list)

twoB = ggplot(tt,aes(x = model, y = x,color = model)) + geom_boxplot(size=1.3)+#geom_line(size=1.3,group = 1)+
  labs(y='Ratio of variance to mean',x="")+theme_bw()+#ylim(0,64)+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single-cell"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  geom_hline(aes(yintercept=(1),linetype = "twodash"),size = 1.5)+
  theme(strip.text.x = element_text(size = 20,face = "bold"),axis.title=element_text(size=20,face = "bold"),
        legend.position = "none",
        #axis.text.x = element_text(size=20,angle=90,face = "bold"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=20,face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  scale_y_continuous(trans="log2")+
  theme(legend.position = "top",
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        plot.margin = unit(c(0,30,0,0), "pt"),
        legend.box.margin=margin(20,20,20,20))+
  guides(color=guide_legend(nrow=3))

ggsave('twoB.png',dpi = 2000,width = 7.53 , height = 5.72)

plt = ggarrange(twoA,twoB,common.legend = TRUE, legend="top")
ggsave('plot/twoAB.png',dpi = 1500,width = 15 , height = 5.72)


meanVar = list()
total = readRDS("./overdispersionParam.rds")
for(i in fourList){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  phi = coef(nls(obs_var ~ obs_mean + obs_mean^2/th, data = dat,start = list(th = 1),trace = F,control = nls.control(warnOnly=T)))
  dat$newVar<- dat$obs_mean+(phi)*(dat$obs_mean)^2

  meanVar[[i]]<-ggplot(dat, mapping = aes(x = log(obs_mean), y = (obs_var)))+geom_point()+
    theme_bw() +labs(x='Gene Mean',y='Variance',title = i)+
    theme(panel.grid =element_blank(),
          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line =element_line(colour = "black"),
          axis.title = element_blank(),
          axis.text = element_text(size = 16,face = "bold"),
          title = element_text(size=18))+
    geom_line(aes(x = log(obs_mean), y = obs_mean), color = "blue",size=1.3)+
    geom_line(aes(x = log(obs_mean), y = newVar), color = "red",size=1.3)+
    scale_y_continuous(trans="log",breaks = trans_breaks("log", function(x) exp(x)),
                       labels = trans_format("log", math_format(.x)))
  rm(tmp,dat)
}
twoC = plot_grid(meanVar[[1]], meanVar[[2]], meanVar[[3]], ncol=3,align = 'v')
twoC = annotate_figure(twoC,left = text_grob("Variance",face = "bold", size = 18,rot =90),
                       bottom = text_grob("Gene Mean",face = "bold", size = 18))
ggsave('twoC.png',dpi = 2000,width = 7.3 , height = 5.72)
png('plot/twoC.png',width = 7.53 , height = 5.72,res=2000)
plot(twoC)
dev.off()


mseVarP = NULL
mseVarNB = NULL
total = readRDS("./overdispersionParam.rds")
for(i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  tmp = total[which(total$model==i),]
  dat$newVar<-dat$obs_mean+(dat$obs_mean)^2*tmp$theta
  mseVarP<-c(mseVarP,mean((log(dat$obs_var) - log(dat$obs_mean))^2))
  mseVarNB<-c(mseVarNB,mean((log(dat$newVar) - log(dat$obs_var))^2))
  
}

mseVar<-data.frame(poisson=mseVarP,NB=mseVarNB)
mseVar<-gather(mseVar,model,Val)
mseVar$model<-factor(mseVar$model,levels = c("poisson","NB"))
mseVar$tech<-factor(rep(list,2),levels = list)
twoD = ggplot(data = mseVar, mapping = aes(x = tech, y = (Val), colour = model,group = model))+ geom_line(size = 1.3)+
  labs(y='MSE of log Variance',x='')+theme_bw()+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  geom_vline(aes(xintercept=c(8.5),linetype = "non-UMI data"),color = "grey",size = 1,key_glyph = "path",show.legend = F)+
  theme(axis.title=element_text(size=20,face = "bold"),
        legend.position = c(0.85, 0.85),
        axis.text.y = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(size = 20,face = "bold",angle = 90,vjust = 0.4),axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=20,face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=17, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))+
  scale_colour_manual(values = c("blue", "red"))+
  scale_y_continuous(trans = 'log',
                     labels = (scientific))


ggsave('./plot/twoD.png',dpi = 2000,width = 7.3 , height = 5.72)

plot_grid(twoA,twoB,twoC,twoD,labels = c('A', 'B','C','D'),nrow=2,align = 'v')


## plot 3 (proportion bar plot)-----
pp = NULL
nb = NULL
zip = NULL
zinb = NULL
for(i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  pp = c(pp,length(which(dat$final_chooseAIC==1))/nrow(dat))
  nb = c(nb,length(which(dat$final_chooseAIC==2))/nrow(dat))
  zip = c(zip,length(which(dat$final_chooseAIC==3))/nrow(dat))
  zinb = c(zinb,length(which(dat$final_chooseAIC==4))/nrow(dat))
}

tmp = data.frame(tech=list, Poisson=pp, NB=nb, ZIP=zip, ZINB=zinb)
tmp$tech = factor(tmp$tech,levels = list)
tmp = gather(tmp,model,val,-tech)
tmp$model = factor(tmp$model,levels = c("Poisson","NB","ZIP","ZINB"))

p1 = ggplot(tmp,aes(x = tech, y = val,fill = model)) + geom_bar(stat="identity",position = position_fill(reverse = TRUE)) + 
  labs(y='Proportion of preferred model',x="")+theme_bw()+#ylim(0,64)+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell data"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  theme(strip.text.x = element_text(size = 18),legend.title= element_blank(),
        legend.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size = 18,angle = 90,hjust=0.5), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'top')+
  # guides(colour = guide_legend(order = 1), 
  #        shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_brewer(palette = "Set2")


### proportion of significant gene
chart<-NULL
for (i in list){
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  overDis<-which(dat$final_PvsNB<(0.05/nrow(dat)))
  zip<-which(dat$PvsZIP<(0.05/nrow(dat)))
  zinb<-which(dat$final_NBvsZINB<(0.05/nrow(dat)))
  zi = which(dat$ZIPvsZINB<(0.05/nrow(dat)))
  
  #tmp = c((overDis),(zip),(zinb),(zinb),(zi))
  tmp<-c(length(overDis),length(zip),length(zinb),length(zi))/nrow(dat)
  chart<-rbind(chart,tmp)
  rm(tmp)
}
chart = as.data.frame(chart)
colnames(chart) = c("PvsNB","PvsZIP","NBvsZINB","ZIPvsZINB")
tmp = gather(chart, test,value)
tmp$model = factor(rep(list,4),levels = list)
tmp$test = factor(tmp$test,levels = c("PvsZIP","PvsNB","NBvsZINB","ZIPvsZINB"))
plt = ggplot(data = tmp, mapping = aes(x = model, y = value, color = test,group = test))+ geom_point(size = 3,stat="identity")+
  geom_line(size=1.3)+
  labs(y='Proportion of significant genes',x='',title = "")+theme_bw()+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell data"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  theme(strip.text.x = element_text(size = 18),legend.title= element_blank(),
        legend.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text.y = element_text(size = 18,angle = 90,hjust=0.5), 
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 18,angle = 90,vjust=0.5),
        plot.title = element_text (hjust = 0.5),title = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'top')+
  # guides(colour = guide_legend(order = 1), 
  #        shape = guide_legend(order = 2))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

plot3 = ggarrange(p1,plt, labels = c("A", "B"),ncol = 1, nrow = 2)

png('plot/plot3.png',width = 30,height = 28,res = 2000,units = 'cm')
plot(plot3)
dev.off()


## new plot 4 -----
res1 = data.frame(matrix(ncol = 5, nrow = 0))
colnames(res1) = c('Poisson','NB','ZIP','ZINB','tech')
for (i in list) {
  print(i)
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  base<-summaryAICNB(dat[which(dat$obs_zerofrac<1),])[1:4]
  rm(dat)
  data<-readRDS(paste0("celltype_AIC/",i,"_celltype.rds"))
  sum = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(sum) = c('Poisson','NB','ZIP','ZINB')
  for(j in 1:length(data)){
    tryCatch({
      dat<-data[[j]]
      tmp = summaryAICNB(dat[which(dat$obs_zerofrac<1),])[1:4]
      print(tmp)
      sum[j,] = tmp/base
      #tmp<-c(tmp,summaryAICNB(dat[which(dat$obs_zerofrac>0),]))
      rm(dat)
    }, error = function(e) {
    })
  }
  sum = apply(sum, 2, function(x) median(x,na.rm=TRUE))
  sum$tech = i
  res1 = rbind(res1,sum)
}
res1[res1==0] = 1
res1 = gather(res1, model, val, -tech)
res1$model = factor(res1$model,levels = c('Poisson','NB','ZIP','ZINB'))
res1$tech = factor(res1$tech,levels=list)
fourA = ggplot(data = res1, mapping = aes(x = tech, y = (val), color = model,group = model))+ geom_point(size = 3,stat="identity")+
  geom_line(size=1.8)+#ylim(-3,3)+
  labs(y='Median ratio on the preferred models \n before and after clustering',x='',title = "")+
  theme_bw()+#ylim(0,64)+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell data"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
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
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(trans = 'log',
                     labels = scales::number_format(accuracy = 0.01L))



res = data.frame(matrix(ncol = 5, nrow = 0))
colnames(res) = c('PvsZIP','PvsNB','NBvsZINB','ZIPvsZINB')
for (i in list) {
  dat<-readRDS(paste0("new2/",i,"_sum.rds"))
  dat<-dat[which(dat$obs_zerofrac<1),]
  overDis<-which(dat$final_PvsNB<(0.05/nrow(dat)))
  zip<-which(dat$PvsZIP<(0.05/nrow(dat)))
  zinb<-which(dat$final_NBvsZINB<(0.05/nrow(dat)))
  zi = which(dat$ZIPvsZINB<(0.05/nrow(dat)))
  base = c(length(zip)/nrow(dat),length(overDis)/nrow(dat),length(zinb)/nrow(dat),length(zi)/nrow(dat))
  rm(dat)
  data<-readRDS(paste0("celltype_AIC/",i,"_celltype.rds"))
  sum = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(sum) = c('PvsZIP','PvsNB','NBvsZINB','ZIPvsZINB')
  for(j in 1:length(data)){
    tryCatch({
      dat<-data[[j]]
      dat<-dat[which(dat$obs_zerofrac<1),]
      overDis<-which(dat$final_PvsNB<(0.05/nrow(dat)))
      zip<-which(dat$PvsZIP<(0.05/nrow(dat)))
      zinb<-which(dat$final_NBvsZINB<(0.05/nrow(dat)))
      zi = which(dat$ZIPvsZINB<(0.05/nrow(dat)))
      tmp = c(length(zip)/nrow(dat),length(overDis)/nrow(dat),length(zinb)/nrow(dat),length(zi)/nrow(dat))
      sum[j,] = tmp/base
      #tmp<-c(tmp,summaryAICNB(dat[which(dat$obs_zerofrac>0),]))
      rm(dat)
    }, error = function(e) {
      
    })
  }
  sum = apply(sum, 2, function(x) median(x,na.rm=TRUE))
  res = rbind(res,sum)
  print(i)
}
res[res == 0]=1 ## just for plotting
res$tech = factor(list,levels = list)
colnames(res) = c('PvsZIP','PvsNB','NBvsZINB','ZIPvsZINB','tech')
res = gather(res, model, val, -tech)
res$model = factor(res$model,levels = c('PvsZIP','PvsNB','NBvsZINB','ZIPvsZINB'))
res$tech = factor(res$tech,levels=list)

fourB = ggplot(data = res, mapping = aes(x = tech, y = (val), color = model,group = model))+ geom_point(size = 3,stat="identity")+
  geom_line(size=1.8)+#ylim(-3,3)+
  labs(y='Median ratio on the significant genes \n before and after clustering',x='',title = "")+
  theme_bw()+#ylim(0,64)+
  geom_vline(aes(xintercept=c(3.5),linetype = "smFISH"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
  geom_vline(aes(xintercept=c(8.5),linetype = "single cell data"),color = "grey",size = 1,key_glyph = "path",show.legend=FALSE)+
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

ggsave("plot/plot4.png",plot = plot4,width = 35,height = 32,
       units = 'cm', dpi = 2000)

