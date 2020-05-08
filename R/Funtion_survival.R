
# Function I. bootstrap for univariate analysis
# 1. statistic parameter used for "boot" function
uni_boot <-function(data,indices,vari){
  d = data[indices,] # d should contain time and status, indices is required for boot function
  print(indices)
  print(length(indices))
  print(length(unique(indices)))
  unifit = coxph(as.formula(paste('Surv(time, status)~', vari))
                ,data = d)
  # return hazard ratio
  summary(unifit)$coef[2]
}

# 2. univariate boot
mul_uni_boot<-function(data,vari,R){
  r1 <- boot(data = data, statistic = uni_boot, R = R,
            vari = vari)
  # "percentile" CI is used; percent for type = "perc"; normal for type="norm"; 
  r2 <- boot.ci(boot.out = r1, type = "perc")$percent[c(4,5)]
  # original value, mean, bias, 95% CI of defined type
  return(c(r1$t0,mean(r1$t),r1$t0-mean(r1$t),r2))
}
# comment: you should always construct a boot statistic firstly before using boot function;
#          for survival analysis, the data you pass to boot should include depedent variables used in statistic parameter (here time and status)
#          there are different types of confidence interval from boot.ci
#          Finally, I would like to recommned using replicate instead of boot function when performing bootstrap

# Function II: Bootstrap used for stepwise
stepwise_boot <-function(mydata,time,status,retu_sta=T){
  # on average, 2/3 samples would be sampled
  ID = sample(1:nrow(mydata),nrow(mydata),replace = T) 
  # train data 
  fit_data = mydata[ID,]
  fit_time = time[ID]
  fit_status = status[ID]
  # out of bag
  tes_data = mydata[-unique(ID),]
  tes_time = time[-unique(ID)]
  tes_status = status[-unique(ID)]
  # train the model
  step_formulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(colnames(fit_data), sep="", collapse="+")))
  modelAll.coxph <- coxph(step_formulas,data = fit_data)
  result.step <- step(modelAll.coxph)
  selected_form <- as.formula(paste('Surv(fit_time,fit_status)~',paste(names(result.step$coefficients),sep="",collapse="+")))
  
  if(retu_sta == T){
    selected.coxph <- coxph(selected_form,data = fit_data)
    risk_rank = predict(selected.coxph,tes_data,type="risk") # we have problem with extreme high or low risk
    #should we delete those points or not? The anser is no since they are part of your result.
    risk_mark = ifelse(risk_rank>median(risk_rank),0,1) # low risk as 1
    # log-rant test statistics
    logran = survdiff(Surv(tes_time,tes_status) ~ risk_mark, rho=0)$chisq
  }else{
    c(names(result.step$coefficients))
  }
}

# Function III: bootstrap for evaluating mature models
model_boot <-function(mydata,time,status,coef.vec,R.r=T){

  # on average, 2/3 samples would be sampled
  ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
  # train data 
  fit_data = mydata[ID,]
  fit_time = time[ID]
  fit_status = status[ID]
  # out of bag
  tes_data = mydata[-unique(ID),]
  tes_time = time[-unique(ID)]
  tes_status = status[-unique(ID)]
  # train the model
  kk = 1
  if(length(coef.vec)>1){ # multiple variables
      formulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(coef.vec, sep="", collapse="+")))
  }else{ # univariate
      print(paste("Variable",coef.vec,"is being processed."))
      formulas <- as.formula(paste('Surv(fit_time,fit_status)~', coef.vec))
  }
  selected.coxph <- coxph(formulas,data = fit_data)
  tt<-tryCatch(coxph(formulas,data = fit_data),error=function(e) e, warning=function(w) w)
  while(is(tt,"warning")){
      print("Not converged")
      # on average, 2/3 samples would be sampled
      ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
      # train data
      fit_data = mydata[ID,]
      fit_time = time[ID]
      fit_status = status[ID]
      # out of bag
      tes_data = mydata[-unique(ID),]
      tes_time = time[-unique(ID)]
      tes_status = status[-unique(ID)]
      # train the model
      selected.coxph <- coxph(formulas,data = fit_data)
      tt<-tryCatch(coxph(formulas,data = fit_data),error=function(e) e, warning=function(w) w)
      kk=kk+1
  }
  print(paste("The process is repeated for",kk,"times."))
  #print(summary(selected.coxph)$coef[,2:3])
  risk_rank = predict(selected.coxph,tes_data,type="risk") # we have problem with extreme high or low risk
  #should we delete those points or not, no, since they are outliers of bootstrap
  risk_mark = ifelse(risk_rank>median(risk_rank),0,1) # low risk as 1
  # log-rant test statistics
  logran = survdiff(Surv(tes_time,tes_status) ~ risk_mark, rho=0)$chisq
  if(R.r==T){
      return(c(logran,summary(selected.coxph)$rsq[1]))
  }else{
      return(logran)
  }
}

# Function IV: univariate analysis
multiple_uni_cox<-function(time,status,All_varData){
    # formulas for each variable
    covariates = colnames(All_varData)
    
    univ_formulas <- sapply(covariates,
    function(x) as.formula(paste('Surv(time, status)~', x)))
    
    univ_models <- lapply(univ_formulas, function(x){coxph(x, data = All_varData)})
    # Extract data
    univ_results <- lapply(univ_models,
    function(x){
        y <- x
        x <- summary(x)
        if(dim(x$coef)[1]==1){                         # we have quantitative variable
            p.value<-signif(x$wald["pvalue"], digits=3)
            wald.test<-signif(x$wald["test"], digits=3)
            beta<-signif(x$coef[1], digits=3);#coeficient beta
            HR <-signif(x$coef[2], digits=3);#exp(beta)
            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
            HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
            logHR_upper <- signif(beta+1.96*x$coef[3],3)
            logHR_lower <- signif(beta-1.96*x$coef[3],3)
            HR1 <- paste0(HR, " (",
            HR.confint.lower, "-", HR.confint.upper, ")")
            # concordance
            cindex <- concordance.index(predict(y),surv.time = time, surv.event = status ,method = "noether")
            res<-c(beta, logHR_lower,logHR_upper,HR1, wald.test, p.value,HR,HR.confint.lower,HR.confint.upper,cindex$c.index,cindex$lower,cindex$upper)
            names(res)<-c("beta", "coef_95%_lower","coef_95%_upper","HR (95% CI for HR)", "wald.test",
            "p.value","HR","HR_lower","HR_upper","Concordance","C_index_lower","C_index_upper")
            #return(exp(cbind(coef(x),confint(x))))
            return(res)
        }else{ # deal with factors
            p.value<-signif(x$wald["pvalue"], digits=3)
            wald.test<-signif(x$wald["test"], digits=3)
            beta<-signif(x$coef[,1], digits=3);#coeficient beta
            HR <-signif(x$coef[,2], digits=3);#exp(beta)
            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
            HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
            logHR_upper <- signif(beta+1.96*x$coef[,3],3)
            logHR_lower <- signif(beta-1.96*x$coef[,3],3)
            # concordance
            cindex <- concordance.index(predict(y),surv.time = time, surv.event = status ,method = "noether")
            HR1 <- paste0(HR, " (",
            HR.confint.lower, "-", HR.confint.upper, ")")
            res<-cbind(beta, logHR_lower,logHR_upper,HR1, wald.test, p.value,HR,HR.confint.lower,HR.confint.upper,cindex$c.index,cindex$lower,cindex$upper)
            colnames(res)<-c("beta", "coef_95%_lower","coef_95%_upper","HR (95% CI for HR)", "wald.test",
            "p.value","HR","HR_lower","HR_upper","Concordance","C_index_lower","C_index_upper")
            return(res)
        }})
    return(univ_results)
}

# Function V: detect outliers && replace outliers
V_outliers <- function(x, na.rm = T, change = F, fold = 1.5){
  qnt <- quantile(x, probs = c(.25,.75),na.rm = na.rm)
  H <- fold*IQR(x, na.rm = na.rm) # the definition of outliers is 1.5*IQR, here I tried fold 1.5, 5, and 10
  if(!change){ # just return summary of outliers
    out_ID = which((x< (qnt[1]-H))|(x >(qnt[2] + H)))
    return(out_ID)
  }else{ # replace outliers with NA
    x[which((x< (qnt[1]-H))|(x >(qnt[2] + H)))] = NA
    x
  }
}

# Function VI: remove correlated variables
RemoveCor<-function(data,C.cri){
  cor.res = signif(cor(data),3)
  abs.cor = abs(cor.res)
  dup.rec = which(abs.cor>C.cri,arr.ind=T)
  
  nco = c(1:dim(data)[2])
  for(i in 1:dim(dup.rec)[1]){
    if(dup.rec[i,1]<=dup.rec[i,2]){
      next
    }    
    if(dup.rec[i,1]%in%nco){
     if(!(dup.rec[i,2]%in%nco)){
       next
     }
     dup.id =match(dup.rec[i,1],nco)
     nco = nco[-dup.id]
   }
  }
  if(length(nco)==1){
    data1=matrix(data[,nco])
    colnames(data1) = colnames(data)[nco]
  }else{
    data1 = data[,nco]
  }
  return(data1)
}

# Function VII: boot all
boot.every <-function(mydata,time,status,myHeat2,Var_RnameV,R.r=F){
  two_values= c()
  # 1. split data: on average, 2/3 samples would be sampled
  ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
  # train data 
  fit_data = mydata[ID,]
  fit_time = time[ID]
  fit_status = status[ID]
  # out of bag
  tes_data = mydata[-unique(ID),]
  tes_time = time[-unique(ID)]
  tes_status = status[-unique(ID)]
  
  # 2. univariate analysis and filtering
  ## full dataset
  univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
  res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
  sig.name = names(res[res<0.01])
  fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
  print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  while(ncol(fit_S1)==0){
    
    # 1. split data: on average, 2/3 samples would be sampled
    ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
    # train data 
    fit_data = mydata[ID,]
    fit_time = time[ID]
    fit_status = status[ID]
    # out of bag
    tes_data = mydata[-unique(ID),]
    tes_time = time[-unique(ID)]
    tes_status = status[-unique(ID)]
    
    # 2. univariate analysis and filtering
    ## full dataset
    univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
    res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
    sig.name = names(res[res<0.01])
    fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
    print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  }  
  
  ## stable dataset
  cov_m = as.vector(apply(myHeat2,1,median)<10)
  stable_list = Var_RnameV[cov_m]
  if(sum(colnames(fit_S1)%in%stable_list)==1){
    SSfit_S1 = as.matrix(fit_S1[,colnames(fit_S1)%in%stable_list])
    colnames(SSfit_S1) = colnames(fit_S1)[which(colnames(fit_S1)%in%stable_list)]
  }else{
    SSfit_S1 = fit_S1[,colnames(fit_S1)%in%stable_list]
  }
  print(paste0(ncol(SSfit_S1)," stable variables with Pval1 < 0.01"))
  while(ncol(SSfit_S1)==0){
  # 1. split data: on average, 2/3 samples would be sampled
  ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
  # train data 
  fit_data = mydata[ID,]
  fit_time = time[ID]
  fit_status = status[ID]
  # out of bag
  tes_data = mydata[-unique(ID),]
  tes_time = time[-unique(ID)]
  tes_status = status[-unique(ID)]
  
  # 2. univariate analysis and filtering
  ## full dataset
  univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
  res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
  sig.name = names(res[res<0.01])
  fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
  print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  
  while(ncol(fit_S1)==0){
    
    # 1. split data: on average, 2/3 samples would be sampled
    ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
    # train data 
    fit_data = mydata[ID,]
    fit_time = time[ID]
    fit_status = status[ID]
    # out of bag
    tes_data = mydata[-unique(ID),]
    tes_time = time[-unique(ID)]
    tes_status = status[-unique(ID)]
    
    # 2. univariate analysis and filtering
    ## full dataset
    univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
    res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
    sig.name = names(res[res<0.01])
    fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
    print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  }  
  
  ## stable dataset
  cov_m = as.vector(apply(myHeat2,1,median)<10)
  stable_list = Var_RnameV[cov_m]
  if(sum(colnames(fit_S1)%in%stable_list)==1){
    SSfit_S1 = as.matrix(fit_S1[,colnames(fit_S1)%in%stable_list])
    colnames(SSfit_S1) = colnames(fit_S1)[which(colnames(fit_S1)%in%stable_list)]
  }else{
    SSfit_S1 = fit_S1[,colnames(fit_S1)%in%stable_list]
  }
  print(paste0(ncol(SSfit_S1)," stable variables with Pval2 < 0.01"))
}

  # 3. filter with correlation analysis
  fit_S2 = RemoveCor(data=cbind(fit_S1,fit_S1),C.cri=0.75)
  coef.vec = colnames(fit_S2)
  print(paste0(ncol(fit_S2)," full variables after correlation"))

  # 4. train the model
  kk = 1
  if(length(coef.vec)>1){ # multiple variables
   formulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(coef.vec, sep="", collapse="+")))
  }else{
    formulas <- as.formula(paste('Surv(fit_time,fit_status)~', coef.vec))
  }
  selected.coxph <- coxph(formulas,data = fit_data)
  tt<-tryCatch(coxph(formulas,data = fit_data),error=function(e) e, warning=function(w) w)


  # 3. filter with correlation analysis
  SSfit_S2 = RemoveCor(data=cbind(SSfit_S1,SSfit_S1),C.cri=0.75)
  SScoef.vec = colnames(SSfit_S2)
  print(paste0(ncol(SSfit_S2)," full variables after correlation"))

  # 4. train the model
  jj = 1
  if(length(SScoef.vec)>1){ # multiple variables
  SSformulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(SScoef.vec, sep="", collapse="+")))
}else{ # univariate
  # print(paste("Variable",fit.coef.vec,"is being processed."))
  SSformulas <- as.formula(paste('Surv(fit_time,fit_status)~', SScoef.vec))
}
  SSselected.coxph <- coxph(SSformulas,data = fit_data)
  ww<-tryCatch(coxph(SSformulas,data = fit_data),error=function(e) e, warning=function(w) w)

  while(is(tt,"warning")|is(ww,"warning")){
  print("Not converged")
  # on average, 2/3 samples would be sampled
  # 1. split data: on average, 2/3 samples would be sampled
  ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
  # train data 
  fit_data = mydata[ID,]
  fit_time = time[ID]
  fit_status = status[ID]
  # out of bag
  tes_data = mydata[-unique(ID),]
  tes_time = time[-unique(ID)]
  tes_status = status[-unique(ID)]
  
  # 2. univariate analysis and filtering
  ## full dataset
  univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
  res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
  sig.name = names(res[res<0.01])
  fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
  print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  
  while(ncol(fit_S1)==0){
    
    # 1. split data: on average, 2/3 samples would be sampled
    ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
    # train data 
    fit_data = mydata[ID,]
    fit_time = time[ID]
    fit_status = status[ID]
    # out of bag
    tes_data = mydata[-unique(ID),]
    tes_time = time[-unique(ID)]
    tes_status = status[-unique(ID)]
    
    # 2. univariate analysis and filtering
    ## full dataset
    univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
    res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
    sig.name = names(res[res<0.01])
    fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
    print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
  }  
  
  ## stable dataset
  cov_m = as.vector(apply(myHeat2,1,median)<10)
  stable_list = Var_RnameV[cov_m]
  if(sum(colnames(fit_S1)%in%stable_list)==1){
    SSfit_S1 = as.matrix(fit_S1[,colnames(fit_S1)%in%stable_list])
    colnames(SSfit_S1) = colnames(fit_S1)[which(colnames(fit_S1)%in%stable_list)]
  }else{
    SSfit_S1 = fit_S1[,colnames(fit_S1)%in%stable_list]
  }
  print(paste0(ncol(SSfit_S1)," stable variables with Pval3 < 0.01"))
  
  while(ncol(SSfit_S1)==0){
    # 1. split data: on average, 2/3 samples would be sampled
    ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
    # train data 
    fit_data = mydata[ID,]
    fit_time = time[ID]
    fit_status = status[ID]
    # out of bag
    tes_data = mydata[-unique(ID),]
    tes_time = time[-unique(ID)]
    tes_status = status[-unique(ID)]
    
    # 2. univariate analysis and filtering
    ## full dataset
    univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
    res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
    sig.name = names(res[res<0.01])
    fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
    print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
    
    while(ncol(fit_S1)==0){
      
      # 1. split data: on average, 2/3 samples would be sampled
      ID = sample(1:nrow(mydata),nrow(mydata),replace = T)
      # train data 
      fit_data = mydata[ID,]
      fit_time = time[ID]
      fit_status = status[ID]
      # out of bag
      tes_data = mydata[-unique(ID),]
      tes_time = time[-unique(ID)]
      tes_status = status[-unique(ID)]
      
      # 2. univariate analysis and filtering
      ## full dataset
      univ_results = multiple_uni_cox(time=fit_time,status=fit_status,All_varData=fit_data)
      res <- sort(t(as.data.frame(univ_results[1:841], check.names = FALSE))[,6])
      sig.name = names(res[res<0.01])
      fit_S1 = fit_data[,match(sig.name,colnames(fit_data))]
      print(paste0(ncol(fit_S1)," variables with Pval < 0.01"))
    }  
    
    ## stable dataset
    cov_m = as.vector(apply(myHeat2,1,median)<10)
    stable_list = Var_RnameV[cov_m]
    if(sum(colnames(fit_S1)%in%stable_list)==1){
      SSfit_S1 = as.matrix(fit_S1[,colnames(fit_S1)%in%stable_list])
      colnames(SSfit_S1) = colnames(fit_S1)[which(colnames(fit_S1)%in%stable_list)]
    }else{
      SSfit_S1 = fit_S1[,colnames(fit_S1)%in%stable_list]
    }
    print(paste0(ncol(SSfit_S1)," stable variables with Pval4 < 0.01"))
  }
  
  # 3. filter with correlation analysis
  fit_S2 = RemoveCor(data=cbind(fit_S1,fit_S1),C.cri=0.75)
  coef.vec = colnames(fit_S2)
  print(paste0(ncol(fit_S2)," full variables after correlation"))
  
  # 4. train the model
  kk = 1
  if(length(coef.vec)>1){ # multiple variables
    formulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(coef.vec, sep="", collapse="+")))
  }else{ # univariate
    # print(paste("Variable",fit.coef.vec,"is being processed."))
    formulas <- as.formula(paste('Surv(fit_time,fit_status)~', coef.vec))
  }
  selected.coxph <- coxph(formulas,data = fit_data)
  tt<-tryCatch(coxph(formulas,data = fit_data),error=function(e) e, warning=function(w) w)
  kk=kk+1
  
  # 3. filter with correlation analysis
  SSfit_S2 = RemoveCor(data=cbind(SSfit_S1,SSfit_S1),C.cri=0.75)
 
  SScoef.vec = colnames(SSfit_S2)
  print(paste0(ncol(SSfit_S2)," stable variables after correlation"))
  
  # 4. train the model
  jj = 1
  if(length(SScoef.vec)>1){ # multiple variables
    SSformulas <- as.formula(paste('Surv(fit_time,fit_status)~', paste(SScoef.vec, sep="", collapse="+")))
  }else{ # univariate
    # print(paste("Variable",fit.coef.vec,"is being processed."))
    SSformulas <- as.formula(paste('Surv(fit_time,fit_status)~', SScoef.vec))
  }
  SSselected.coxph <- coxph(SSformulas,data = fit_data)
  ww<-tryCatch(coxph(SSformulas,data = fit_data),error=function(e) e, warning=function(w) w)
  jj = jj+1
}
  
  print(paste("The process is repeated for",kk,"times and",jj,"times"))
  #print(summary(selected.coxph)$coef[,2:3])
  risk_rank = predict(selected.coxph,tes_data,type="risk") # we have problem with extreme high or low risk
  #should we delete those points or not, no, since they are outliers of bootstrap
  risk_mark = ifelse(risk_rank>median(risk_rank),0,1) # low risk as 1
  # log-rant test statistics
  logran = survdiff(Surv(tes_time,tes_status) ~ risk_mark, rho=0)$chisq
  two_values[1] = logran
  
  SSrisk_rank = predict(SSselected.coxph,tes_data,type="risk") # we have problem with extreme high or low risk
  #should we delete those points or not, no, since they are outliers of bootstrap
  SSrisk_mark = ifelse(SSrisk_rank>median(SSrisk_rank),0,1) # low risk as 1
  # log-rant test statistics
  SSlogran = survdiff(Surv(tes_time,tes_status) ~ SSrisk_mark, rho=0)$chisq
  two_values[2] = SSlogran
  print(paste0(two_values[1],"and",two_values[2]))
  return(two_values)
}
  
  
# Function VIII: univariate analysis for classification
multiple_uni_anova <-function(featureM,Target){
    sumcr.mva1 <- manova(as.matrix(featureM) ~ factor(Target))
    sum.mano = summary.aov(sumcr.mva1)
    pval.mano = sapply(sum.mano,function(x) x$"Pr(>F)"[1]) # pvalue for each variable
    names(pval.mano) = colnames(featureM)
    return(pval.mano)
}



 
  
  
  
  
  


















