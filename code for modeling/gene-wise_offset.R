
#source("function.R")
ofset<-colSums(dat)
sum_stat <- data.frame(rownames(dat))

# calculate overall statistics
sum_stat$obs_mean<-NA
sum_stat$obs_var<-NA
sum_stat$obs_zerofrac<-NA

sum_stat[,paste0(c('P','NB','ZIP','ZINB'),'_aic')]<-NA
sum_stat[,paste0(c('P','NB','ZIP','ZINB'),'_de')]<-NA
sum_stat[,paste0(c('P','NB','ZIP','ZINB'),'_exp_zerofrac')]<-NA
sum_stat[,c('chooseAIC','chooseDe')]<-NA


## get poisson related statistics
#reference:https://data.princeton.edu/wws509/notes/countmoments

number_of_cores <- parallel::detectCores() - 1
#number_of_cores <- 5
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)

stat<-foreach(i = 1:nrow(dat), .combine = rbind,.packages = 'Matrix')  %dopar%   stat_sum(dat[i,])
P_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS'))  %dopar%   po_fun(dat[i,],ofset)
print("p finish")
NB_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','MASS'))  %dopar%   nb_fun(dat[i,],ofset)
print("NB finish")
ZIP_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zip_fun(dat[i,],ofset)
print("ZIP finish")
ZINB_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zinb_fun(dat[i,],ofset)
print("ZINB finish")
# ZIPcov_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zip_fun_cov(dat[i,],ofset)
# ZINBcov_info<-foreach(i=1:nrow(dat), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zinb_fun_cov(dat[i,],ofset)
testPValue<-foreach(i=1:nrow(dat),.combine = rbind,.packages = c('MASS','lmtest','Matrix','pscl','vcdExtra')) %dopar% test(dat[i,],ofset)

parallel::stopCluster(clusters)

sum_stat$obs_mean<-stat$mean
sum_stat$obs_var<-stat$var
sum_stat$obs_zerofrac<-stat$zerofrac

sum_stat$P_aic <- P_info$aic
sum_stat$P_de <- P_info$de
sum_stat$P_exp_zerofrac <- exp(-P_info$mu)
sum_stat$P_mean <- P_info$mu
sum_stat$P_var <- P_info$mu

sum_stat$NB_aic <- NB_info$aic
sum_stat$NB_de <- NB_info$de
sum_stat$NB_exp_zerofrac <- (1+NB_info$alpha*NB_info$mu)^(-1/NB_info$alpha)
sum_stat$NB_mean <- NB_info$mu
sum_stat$NB_var <- NB_info$mu*(1+NB_info$alpha*NB_info$mu)

sum_stat$ZIP_aic <- ZIP_info$aic
sum_stat$ZIP_de <- ZIP_info$de
sum_stat$ZIP_exp_zerofrac <- ZIP_info$pi+(1-ZIP_info$pi)*exp(-ZIP_info$mu)
sum_stat$ZIP_mean <- (1-ZIP_info$pi)*ZIP_info$mu
sum_stat$ZIP_var <- (1-ZIP_info$pi)*ZIP_info$mu*(1+ZIP_info$mu*ZIP_info$pi)

sum_stat$ZINB_aic <- ZINB_info$aic
sum_stat$ZINB_de <- ZINB_info$de
sum_stat$ZINB_exp_zerofrac <- ZINB_info$pi+(1-ZINB_info$pi)*(1+ZINB_info$alpha*ZINB_info$mu)^(-1/ZINB_info$alpha)
sum_stat$ZINB_mean <- (1-ZINB_info$pi)*ZINB_info$mu
sum_stat$ZINB_var <- (1-ZINB_info$pi)*ZINB_info$mu*(1+ZINB_info$mu*(ZINB_info$pi+ZINB_info$alpha))

sum_stat$chooseAIC<-apply(sum_stat[,paste0(c('P','NB','ZIP','ZINB'),'_aic')],1,chooseAIC)
sum_stat$chooseDe<-apply(sum_stat[,paste0(c('P','NB','ZIP','ZINB'),'_de')],1,chooseDe)

sum_stat$PvsNB<-testPValue$p1
sum_stat$inflation<-testPValue$p2
sum_stat$PvsZIP<-testPValue$p3
sum_stat$NBvsZINB<-testPValue$p4
sum_stat$ZIPvsZINB<-testPValue$p5

sum_stat<-sum_stat %>% remove_rownames %>% column_to_rownames(var=names(sum_stat)[1])

