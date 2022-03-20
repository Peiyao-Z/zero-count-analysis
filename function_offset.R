stat_sum<-function(arr){
  arr = as.numeric(arr)
  return(data.frame(mean = mean(arr),var = var(arr),zerofrac = sum(arr==0)/length(arr)))
}

po_fun<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    po = glm(arr ~ offset(log(ofset)),family = "poisson")
    #po = glm(arr ~ 1,family = "poisson")
    de_po = sum(resid(po)^2)/po$df.residual
    aic_po = AIC(po)
    mu_po = mean(exp(log(ofset)+rep(coef(po),length(ofset))))
    return(data.frame(aic = aic_po,de = de_po,mu = mu_po))
  },error=function(e){
    return(data.frame(aic = NA,de = NA,mu = NA))
  })
}


nb_fun<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    nb = glm.nb(arr ~ 1+offset(log(ofset)))
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

zip_fun<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    zpo = zeroinfl(arr ~ 1+offset(log(ofset))|1,dist = "poisson")
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

# zip_fun_cov<-function(arr,ofset){
#   arr = as.numeric(arr)
#   tryCatch({
#     #ofset = rep(sum(arr),length(arr))
#     zpo = zeroinfl(arr ~ 1+offset(log(ofset))|1+log(ofset),dist = "poisson")
#     #zpo = zeroinfl(arr ~ 1|1,offset = log(ofset),dist = "poisson")
#     de_zpo = sum(resid(zpo)^2)/zpo$df.residual
#     aic_zpo = AIC(zpo)
#     #mu_zpo = exp(coef(zpo)[1])
#     mu_zpo = mean(exp(log(ofset)+rep(coef(zpo)[1],length(ofset))))
#     pi_zpo = exp(coef(zpo)[3])/(1+exp(coef(zpo)[3]))
#     return(data.frame(aic = aic_zpo,de = de_zpo,mu = mu_zpo,pi = pi_zpo))
#   },error=function(e){
#     return(data.frame(aic = NA,de = NA,mu = NA,pi = NA))
#   })
# }

zinb_fun<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    znb = zeroinfl(arr ~ 1+offset(log(ofset))|1,dist = "negbin")
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

# zinb_fun_cov<-function(arr,ofset){
#   arr = as.numeric(arr)
#   tryCatch({
#     #ofset = rep(sum(arr),length(arr))
#     znb = zeroinfl(arr ~ 1+offset(log(ofset))|1+log(ofset),dist = "negbin")
#     de_znb = sum(resid(znb)^2)/znb$df.residual
#     aic_znb = AIC(znb)
#     mu_znb = mean(exp(log(ofset)+rep(coef(znb)[1],length(ofset))))
#     pi_znb = exp(coef(znb)[3])/(1+exp(coef(znb)[3]))
#     alpha_znb = 1/znb$theta
#     return(data.frame(aic = aic_znb,de = de_znb,mu = mu_znb,pi = pi_znb, alpha = alpha_znb))
#   },error=function(e){
#     return(data.frame(aic = NA,de = NA,mu = NA,pi = NA,alpha=NA))
#   })
# }


chooseAIC<-function(arr){
  if(sum(!is.na(arr))==0){
    return(NA)
  }
  else{
    t = as.numeric(which.min(arr))
    return(t)
  }
}
chooseDe<-function(arr){
  if(sum(!is.na(arr))==0){
    return(NA)
  }
  else{
    t = as.numeric(which.min(arr-1))
    return(t)
  }
}

test<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
    #ofset = rep(sum(arr),length(arr))
    po = glm(arr ~ offset(log(ofset)),family = "poisson")
    nb = glm.nb(arr ~ offset(log(ofset)))
    zip = zeroinfl(arr ~ 1|1,offset = log(ofset),dist = "poisson")
    zinb = zeroinfl(arr ~ 1|1,offset = log(ofset),dist = "negbin")
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

vuong<-function(arr,ofset){
  arr = as.numeric(arr)
  tryCatch({
  po = glm(arr ~ offset(log(ofset)),family = "poisson")
  #nb = glm.nb(arr ~ offset(log(ofset)))
  zip1 = zeroinfl(arr ~ 1|1,offset = log(ofset),dist = "poisson")
  zip2 = zeroinfl(arr ~ 1|1+log(ofset),offset = log(ofset),dist = "poisson")
  #zinb = zeroinfl(arr ~ 1|1+log(ofset),offset = log(ofset),dist = "negbin")
  tt = vuongtest(po,zip2)
  return(data.frame(A = tt$p_LRT$A, B = tt$p_LRT$B))
  #vuong(nb,zinb)
  },error=function(e){
    return(data.frame(A = NA, B = NA))
  })
}
