po_fun<-function(arr,cov,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    po = glm(arr ~ cov + offset(log(ofset)),family = "poisson")
    de_po = sum(resid(po)^2)/po$df.residual
    aic_po = AIC(po)
    mu_po = mean(exp(log(ofset)+rep(coef(po),length(ofset))))
    return(data.frame(aic = aic_po,de = de_po,mu = mu_po))
  },error=function(e){
    return(data.frame(aic = NA,de = NA,mu = NA))
  })
}

nb_fun<-function(arr,cov,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    nb = glm.nb(arr ~ cov+offset(log(ofset)))
    de_nb = sum(resid(nb)^2)/nb$df.residual
    aic_nb = AIC(nb)
    #mu_nb = exp(coef(nb))
    mu_nb = mean(exp(log(ofset)+rep(coef(nb),length(ofset))))
    alpha_nb = 1/nb$theta
    return(data.frame(aic = aic_nb,de = de_nb,mu = mu_nb,alpha = alpha_nb))
  },error=function(e){
    return(data.frame(aic = NA,de = NA,mu = NA,alpha = NA))
  })
}

zip_fun<-function(arr,cov,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    zpo = zeroinfl(arr ~ cov+offset(log(ofset))|1,dist = "poisson")
    de_zpo = sum(resid(zpo)^2)/zpo$df.residual
    aic_zpo = AIC(zpo)
    #mu_zpo = exp(coef(zpo)[1])
    mu_zpo = mean(exp(log(ofset)+rep(coef(zpo)[1],length(ofset))))
    pi_zpo = exp(coef(zpo)[2])/(1+exp(coef(zpo)[2]))
    return(data.frame(aic = aic_zpo,de = de_zpo,mu = mu_zpo,pi = pi_zpo))
  },error=function(e){
    return(data.frame(aic = NA,de = NA,mu = NA,pi = NA))
  })
}

zinb_fun<-function(arr,cov,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    znb = zeroinfl(arr ~ cov+offset(log(ofset))|1,dist = "negbin")
    de_znb = sum(resid(znb)^2)/znb$df.residual
    aic_znb = AIC(znb)
    mu_znb = mean(exp(log(ofset)+rep(coef(znb)[1],length(ofset))))
    pi_znb = exp(coef(znb)[2])/(1+exp(coef(znb)[2]))
    alpha_znb = 1/znb$theta
    return(data.frame(aic = aic_znb,de = de_znb,mu = mu_znb,pi = pi_znb, alpha = alpha_znb))
  },error=function(e){
    return(data.frame(aic = NA,de = NA,mu = NA,pi = NA,alpha=NA))
  })
}

chooseAIC<-function(arr){
  if(sum(!is.na(arr))==0){
    return(NA)
  }
  else{
    t = as.numeric(which.min(arr))
    return(t)
  }
}

test<-function(arr, cov, ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    po = glm(arr ~ cov+offset(log(ofset)),family = "poisson")
    nb = glm.nb(arr ~ cov+offset(log(ofset)))
    zip = zeroinfl(arr ~ cov+offset(log(ofset))|1,dist = "poisson")
    zinb = zeroinfl(arr ~ cov+offset(log(ofset))|1,dist = "negbin")
    p1 = lrtest(po,nb)$`Pr(>Chisq)`[2]
    p2 = zero.test(as.numeric(arr))$prob
    p3 = lrtest(po,zip)$`Pr(>Chisq)`[2]
    p4 = lrtest(nb,zinb)$`Pr(>Chisq)`[2]
    p5 = lrtest(zip,zinb)$`Pr(>Chisq)`[2]
    return(data.frame(p1 = p1, p2 = p2, p3 = p3,p4 = p4,p5 = p5))
  },error=function(e){
    return(data.frame(p1 = NA, p2 = NA, p3 = NA, p4 = NA,p5 = NA))
  })
}


#source("function.R")
ofset<-colSums(spatial_count)
sum_stat <- data.frame(rownames(spatial_count))

# calculate overall statistics
sum_stat[,c('P_aic')]<-NA

number_of_cores <- parallel::detectCores() - 1
# number_of_cores <- 20
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)

P_info<-foreach(i=1:nrow(spatial_count), .combine = rbind,.packages = c('Matrix','MASS'))  %dopar%   po_fun(spatial_count[i,], cov ,ofset)
print("p finish")
NB_info<-foreach(i=1:nrow(spatial_count), .combine = rbind,.packages = c('Matrix','MASS'))  %dopar%   nb_fun(spatial_count[i,], cov ,ofset)
print("NB finish")
ZIP_info<-foreach(i=1:nrow(spatial_count), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zip_fun(spatial_count[i,], cov ,ofset)
print("ZIP finish")
ZINB_info<-foreach(i=1:nrow(spatial_count), .combine = rbind,.packages = c('Matrix','pscl'))  %dopar%   zinb_fun(spatial_count[i,], cov ,ofset)
print("ZINB finish")
testPValue<-foreach(i=1:nrow(spatial_count),.combine = rbind,.packages = c('MASS','lmtest','Matrix','pscl','vcdExtra')) %dopar% test(spatial_count[i,], cov ,ofset)

parallel::stopCluster(clusters)

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

sum_stat$PvsNB<-testPValue$p1
sum_stat$inflation<-testPValue$p2
sum_stat$PvsZIP<-testPValue$p3
sum_stat$NBvsZINB<-testPValue$p4
sum_stat$ZIPvsZINB<-testPValue$p5

