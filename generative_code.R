# This file contains all of the functions needed for reproducing the results of the paper...

library(survival)
library(riskRegression)
library(prodlim)
library(Matrix)
library(matrixStats)
library(glmnet)
library(multcomp)
library(dplyr)
library(stringr)
library(rlang)
library(entropy)

## Rank Entropy - the uncertainty at each rank as an alternative to the rank variance for each exposure.

rankEntropy <- function(d, ranks, P = 1000, epsilon = 0){ # Shannon Entropy for the likelihood to be selected if we grab one observation having rank<=d at random. Is exactly designed for the 'multiple bin case' with counts.
  
  nL <- nrow(ranks)
  nVar <- ncol(ranks)
  
  Entropy <- matrix(nrow = P, ncol = d)
  
  for (p in 1:P){
    ranksR <- residual.rank.randomizer(ranks)
    for (i in 1:d){
      SPcounts <- colSums(ranksR==i)
      
      # epsilon adjustment, when epsilon = 0, the below calculation is not affected.
      SPcounts <- SPcounts[(SPcounts / sum(SPcounts)) > epsilon] 
      
      # entropy calculation
      Entropy[p,i] <- entropy(y = SPcounts)
    }
  }
  
  EntropyList <- colMeans(Entropy)
  return(EntropyList)
}

EntropySum <- function(SP, nL, nVar){
  
  Xi <- 0:nL
  Inners <- numeric(0)
  
  for (i in 1:length(SP)){
    Inners <- c(Inners,binomCoef(nL,Xi)*SP[i]^(Xi)*(1-SP[i])^(nL-Xi)*cumsumLog(nL)) 
  }
  
  fS <- sum(Inners)
  return(fS)
}

binomCoef <- Vectorize(function(n,p){
  if (n == p | p == 0){
    return(1)
  }
  if (p>n){
    return(0)
  }
  else{
    
    nM <- sum(log(1:n))
    dM <- sum(log(1:(n-p))) + sum(log(1:p))
    logFac <- nM - dM
    
    return(exp(logFac))
  }
})

cumsumLog <- function(L){ # L is the upper integer
  
  rL <- c(0)
  for (i in 1:L){
    rL <- c(rL, sum(log(1:i)))
  }
  
  return(rL)
}

# Rank assignment with Wald test

rank.assignment <- function(pollutants, SHbeta){ # input: pollutant names and list of model outputs
  
  ordersSH <- numeric(0) # KEEP FORMAT, INTERACTIONS NOW REPRESENT THEIR 3-EXPOSURE-FAMILIES. MAINS REPRESENT THEIR 1-EXPOSURE-FAMILIES.
  pair.list <- pairwise.name.generator(pollutants)
  
  for (i in (1:length(SHbeta))){
    
    beta <- SHbeta[[i]]$coef #coef(SHbeta[[i]])
    beta <- beta[names(beta) %in% pair.list]
    
    VBeta <- SHbeta[[i]]$vcov #vcov(SHbeta[[i]])
    VBeta <- VBeta[rownames(VBeta) %in% pair.list, colnames(VBeta) %in% pair.list]

    Fam.pvals <- rep(NA,length(names(beta)))
    names(Fam.pvals) <- names(beta)
    
    selected.mains <- names(Fam.pvals)[names(Fam.pvals) %in% pollutants]
    for (j in selected.mains){
      
      k <- which(names(Fam.pvals)==j)
      
      R <- rep(0,length(beta))
      R[k] <- 1
      Q = 1
      
      tSize <- (R%*%beta)^2 * solve(R%*%(VBeta)%*%R)
      Fam.pvals[k] <- pchisq(tSize,df=Q,lower.tail = FALSE)
      
    }
    
    pairs.in.list <- names(Fam.pvals)[!names(Fam.pvals) %in% pollutants]
    
    for(j in pairs.in.list){
      
      k <- which(names(Fam.pvals)==j)
      
      loc <- gregexpr(":",names(Fam.pvals)[k])[[1]][1]
      c1 <- substr(names(Fam.pvals)[k],1,loc-1)
      c2 <- substr(names(Fam.pvals)[k],loc+1,nchar(names(Fam.pvals)[k]))
      
      famIndices <- c(which(names(Fam.pvals) %in% c(c1,c2)), k)
      
      # Describing the hypothesis
      R <- matrix(rep(0,length(beta)),ncol = length(beta), nrow = length(famIndices))
      for (r in 1:length(famIndices)){
        R[r,famIndices[r]] <- 1
      }
      # Degrees of freedom
      Q = nrow(R)
      
      # Test statistic and p value
      tSize <- t((R%*%beta)) %*% solve(R%*%VBeta%*%t(R)) %*% (R%*%beta)
      Fam.pvals[k] <- pchisq(tSize, df = Q, lower.tail = FALSE) # df = length(famIndices)
      
    }
    
    # ... Now to fill in ranks to pollutants, with missings where we have not ranked a pollutant...
    # Following the original idea, now interactions represent the 3 member families, and the main effects represent 1 member families.
    
    full.pval.vec <- rep(1,length(pair.list))
    names(full.pval.vec) <- pair.list
    
    full.pval.vec[pair.list %in% names(Fam.pvals)] <- Fam.pvals
    
    ord <- rank(full.pval.vec)
    ord[!names(ord) %in% names(Fam.pvals)] <- NA # Either not in contention or not selected when in contention!
    ordersSH <- rbind(ordersSH,ord)
    
  }
  
  return(ordersSH)
  
}


## Interactive Lasso

# Helper functions

pairwise.name.generator <- function(string.list){
  form <- numeric(0)
  
  for (i in 1:length(string.list)){
    form <- c(form,string.list[i])
  }
  
  for (i in 1:length(string.list)){
    for (j in (1:length(string.list))[-(1:i)]){
      form <- c(form,str_c(string.list[i],":",string.list[j]))
    }
  }
  return(form)
}

pairwise.paster <- function(string.list){
  form <- numeric(0)
  
  for (i in 1:length(string.list)){
    form <- str_c(form,string.list[i],"+")
  }
  
  for (i in 1:length(string.list)){
    for (j in (1:length(string.list))[-(1:i)]){
      form <- str_c(form,string.list[i],"*",string.list[j],"+")
    }
  }
  form <- substring(form,1,nchar(form)-1)
  return(form)
}

main.paster <- function(string.list){
  form <- numeric(0)
  
  for (i in 1:length(string.list)){
    form <- str_c(form,string.list[i],"+")
  }
  form <- substring(form,1,nchar(form)-1)
  return(form)
}

string.contains <- function(s,l){
  ind.matrix <- matrix(rep(0,length(s)*length(l)),nrow=length(s))
  for(i in 1:length(s)){
    for(j in 1:length(l)){
      ind.matrix[i,j] <- grepl(s[i],l[j])
    }
  }
  return(ind.matrix)
}

which.string.contains <- function(s,l){
  ind <- numeric(0)
  for(i in 1:length(s)){
    ind[i] <- which(l==s[i])
  }
  return(ind)
}


# Function itself

iLassoBoot <- function(jsurvdata, socios, pollutants, lambda.list = NULL, event = "ASTHMA_H"){
  
  datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
  name.select <- pairwise.name.generator(pollutants)
  pair.list <- pairwise.paster(pollutants)
  
  sampled_ids <- sample(1:nrow(datamatrix),replace=T)
  boot_data <- datamatrix[sampled_ids,]
  boot_data[,pollutants] <- t((t(boot_data[,pollutants] ) - colMeans(boot_data[,pollutants] ))/colSds(boot_data[,pollutants] )) # standardizing a priori, as to not standardize the interaction terms themselves (thereby introducing varying scales), but only the separate terms.
  jsurvboot <- jsurvdata[sampled_ids,]
  jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants
  
  strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))
  if (event == "ASTHMA_H"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)
  }
  
  if (event == "death"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$death)),strata = strata.indicator)
  }
  
  # Two step hierarchical lasso with pairwise interactions (including only the main effects lasso as well)
  
  # May want to split the data into two parts to use for the separate steps. Avoids potential overfitting due to cross-validation on the same data twice.
  
  # For using the appropriate weights (adaptive lasso)
  
  lasso_model_select <- glmnet(x=boot_data,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = c(rep(0,length(socios)),rep(1,length(pollutants))))
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  
  minLambda <- lasso_model_select$lambda[which.min(AIC)]
  beta_select <- lasso_model_select$beta[,which.min(AIC)]
  beta_select_main <- beta_select[names(beta_select) %in% pollutants]
  beta.indices <- (beta_select_main != 0)
  
  if (sum(beta_select_main)>0){
    
    Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants),"+",coalesce(pairwise.paster(pollutants[beta.indices]),"1")),sep=""))
    MM <- cbind(boot_data[,1:length(socios)],model.matrix(Pformula,jsurvboot)[,-1])
    
    # We now ensure that the previously selected main effects are kept and we include only potential pairwise interactions associated with previous main effects, but we allow for more main effects being selected still.
    
    pen.index <- c(rep(0,length(socios)),c(1-as.numeric(beta.indices),rep(1,sum(1:(sum(beta.indices)-1))-(sum(beta.indices)==1)))) # Removing the 1 for no. interactions that happens when exactly one pollutant is chosen!
    
    lasso_model_select <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = pen.index)
    
    AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
    minLambda <- lasso_model_select$lambda[which.min(AIC)]
    beta_select <- lasso_model_select$beta[,which.min(AIC)]
    
  }
  
  ## For extracting wald tests: We take the coxph output using the selected beta
  lp_exp <- paste(names(beta_select)[beta_select!=0],collapse = "+")
  if (lp_exp == ""){
    lp_exp <- 1
  }
  Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")",paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
  Routput <- coxph(Rformula,data=jsurvboot, x = FALSE, id = 1:nrow(jsurvboot), y = FALSE)
  
  Foutput <- list("coef" = coef(Routput), "vcov" = vcov(Routput))
  
  return(list(Foutput)) 
  
}

nonAssociationBootstrap <- function(jsurvdata, socios, pollutants, lambda.list = NULL, event = "ASTHMA_H", P = 200){
  
  datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
  name.select <- pairwise.name.generator(pollutants)
  pair.list <- pairwise.paster(pollutants)
  
  sampled_ids <- sample(1:nrow(datamatrix),replace=T)
  boot_data <- datamatrix[sampled_ids,]
  boot_data[,pollutants] <- t((t(boot_data[,pollutants] ) - colMeans(boot_data[,pollutants] ))/colSds(boot_data[,pollutants] )) # standardizing a priori, as to not standardize the interaction terms themselves (thereby introducing varying scales), but only the separate terms.
  jsurvboot <- jsurvdata[sampled_ids,]
  jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants
  
  strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))
  if (event == "ASTHMA_H"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)
  }
  
  if (event == "death"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$death)),strata = strata.indicator)
  }
  
  # For using the appropriate weights (adaptive lasso)
  
  pen.index.factor <- 1
  
  # New stuff with AIC instead:
  lasso_model_select <- glmnet(x=boot_data,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = c(rep(0,length(socios)),rep(1,length(pollutants))*pen.index.factor))
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  
  minLambdaMain <- lasso_model_select$lambda[which.min(AIC)]
  beta_select <- lasso_model_select$beta[,which.min(AIC)]
  beta_select_main <- beta_select[names(beta_select) %in% pollutants]
  beta.indices <- (beta_select_main != 0)
  
  Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants),"+",coalesce(pairwise.paster(pollutants[beta.indices]),"1")),sep=""))
  MM <- cbind(boot_data[,1:length(socios)],model.matrix(Pformula,jsurvboot)[,-1])
  
  # We now ensure that the previously selected main effects are kept and we include only potential pairwise interactions associated with previous main effects, but we allow for more main effects being selected still.
  pen.index <- c(rep(0,length(socios)),c(1-as.numeric(beta.indices),rep(1,sum(1:(sum(beta.indices)-1))-(sum(beta.indices)==1))))
  
  lasso_model_select <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = pen.index)
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  minLambdaInt <- lasso_model_select$lambda[which.min(AIC)]
  
  # Bootstrap ranked lists when no association
  
  beta_select_NA <- matrix(nrow=P, ncol = length(name.select))
  
  Rout <- list()
  p <- 1
  while (p <= P){
    
    NA_lasso <- try(single.Perm.SH(boot_data = boot_data, jsurvStrat = jsurvStrat, minLambdaMain = minLambdaMain, minLambdaInt = minLambdaInt, pollutants = pollutants, socios = socios, jsurvboot = jsurvboot, event = event), silent = TRUE)
    
    if(inherits(NA_lasso, "try-error")){
      next 
    }
    
    Rout[[p]] <- NA_lasso
    
    p <- p + 1
    
  }
  
  return(Rout) 
  
}

single.Perm <- function(jsurvStrat, MM, pen.index, minLambda){
  sample.indices <- sample(1:nrow(jsurvStrat),replace=FALSE)
  
  NAstrat <- jsurvStrat[sample.indices,]
  
  NA_lasso <- glmnet(x=MM,y=NAstrat,family="cox",standardize=FALSE, penalty.factor = pen.index, lambda = minLambda)
  
  return(NA_lasso)
}


single.Perm.SH <- function(boot_data,jsurvStrat, minLambdaMain, minLambdaInt, pollutants, socios, jsurvboot, event){
  
  sample.indices <- sample(1:nrow(jsurvStrat),replace=FALSE)
  
  NAstrat <- jsurvStrat[sample.indices,]
  
  jsurvbootNA <- jsurvboot
  jsurvbootNA$age_end_H <- jsurvboot$age_end_H[sample.indices]
  jsurvbootNA$SEX <- jsurvboot$SEX[sample.indices]
  jsurvbootNA$C_KOM <- jsurvboot$C_KOM[sample.indices]
  jsurvbootNA[,event] <- jsurvboot[sample.indices,event]
  
  PformulaMain <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants)),sep=""))
  MMmain <- cbind(boot_data[,1:length(socios)],model.matrix(PformulaMain,jsurvboot)[,-1])
  pen.index <- c(rep(0,length(socios)), rep(1,length(pollutants)))
  NA_lasso_Main <- glmnet(x=MMmain,y=NAstrat,family="cox",standardize=FALSE, penalty.factor = pen.index, lambda = minLambdaMain)
  
  beta_select <- NA_lasso_Main$beta[,1]
  beta_select_main <- beta_select[names(beta_select) %in% pollutants]
  beta.indices <- (beta_select_main != 0)
  
  if (sum(beta_select_main)>0){
    
    PformulaInt <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants),"+",coalesce(pairwise.paster(pollutants[beta.indices]),"1")),sep=""))
    MMint <- cbind(boot_data[,1:length(socios)],model.matrix(PformulaInt,jsurvboot)[,-1])
    pen.index <- c(rep(0,length(socios)),c(1-as.numeric(beta.indices),rep(1,sum(1:(sum(beta.indices)-1))-(sum(beta.indices)==1)))) 
    NA_lasso_int <- glmnet(x=MMint,y=NAstrat,family="cox",standardize=FALSE, penalty.factor = pen.index, lambda = minLambdaInt)
    
    ## For extracting coxph object:
    lp_exp <- paste(rownames(coef(NA_lasso_int))[coef(NA_lasso_int)[,1]!=0],collapse = "+", sep = "")
    if (lp_exp == ""){
      lp_exp <- 1
    }
    Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")", paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
    Routput <- coxph(Rformula,data=jsurvbootNA, x = FALSE, id = 1:nrow(jsurvbootNA), y = FALSE)
    
  }
  
  if (sum(beta_select_main)==0){
    
    lp_exp <- paste(rownames(coef(NA_lasso_Main))[coef(NA_lasso_Main)[,1]!=0],collapse = "+", sep = "")
    if (lp_exp == ""){
      lp_exp <- 1
    }
    Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")", paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
    Routput <- coxph(Rformula,data=jsurvbootNA, x = FALSE, id = 1:nrow(jsurvbootNA), y = FALSE)
    
    
  }
  
  Foutput <- list("coef" = coef(Routput), "vcov" = vcov(Routput))
  
  return(Foutput)
}


## Shadow functions:

iLassoBootShadow <- function(jsurvdata, socios, pollutants, event = "ASTHMA_H"){
  
  datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
  name.select <- pairwise.name.generator(pollutants)
  
  sampled_ids <- sample(1:nrow(datamatrix),replace=TRUE) # bootstrap
  
  boot_data <- datamatrix[sampled_ids,]
  shadowIds <- sample(1:nrow(boot_data),replace=FALSE) # permutation
  
  shadowPollutants <- boot_data[shadowIds,pollutants]
  shadows <- str_c("Shadow", 1:length(pollutants))
  colnames(shadowPollutants) <- shadows
  
  shadow.select <- pairwise.name.generator(shadows)
  
  boot_data <- cbind(boot_data, shadowPollutants)
  
  boot_data[,c(pollutants,shadows)] <- t((t(boot_data[,c(pollutants,shadows)] ) - colMeans(boot_data[,c(pollutants,shadows)]))/
                                           colSds(boot_data[,c(pollutants,shadows)] )) 
  
  jsurvboot <- jsurvdata[sampled_ids,]
  jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants
  jsurvboot <- cbind(jsurvboot, shadowPollutants)
  
  strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))
  if (event == "ASTHMA_H"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)
  }
  
  if (event == "death"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$death)),strata = strata.indicator)
  }
  
  # Two step hierarchical lasso with pairwise interactions (including only the main effects lasso as well)
  
  lasso_model_select <- glmnet(x=boot_data,y=jsurvStrat,family="cox",standardize=FALSE, 
                               penalty.factor = c(rep(0,length(socios)),rep(1,length(c(pollutants,shadows)))))
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  
  minLambda <- lasso_model_select$lambda[which.min(AIC)]
  beta_select <- lasso_model_select$beta[,which.min(AIC)]
  beta_select_main_p <- beta_select[names(beta_select) %in% c(pollutants)]
  beta_select_main_s <- beta_select[names(beta_select) %in% c(shadows)]
  beta.indices.p <- (beta_select_main_p != 0)
  beta.indices.s <- (beta_select_main_s != 0)
  
  
  if (sum(c(beta_select_main_p,beta_select_main_s))>0){
    
    # SH
    
    Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants), "+",coalesce(pairwise.paster(pollutants[beta.indices.p]),"1"),"+", main.paster(shadows),"+",coalesce(pairwise.paster(shadows[beta.indices.s]),"1")),sep=""))
    MM_pol <- model.matrix(Pformula,jsurvboot)[,-1]
    MM_pol <- MM_pol[,order(match(colnames(MM_pol),c(name.select,shadow.select)))]
    MM <- cbind(boot_data[,1:length(socios)],MM_pol)
    
    pen.index.ps <- c(1-as.numeric(beta.indices.p), rep(1,sum(1:(sum(beta.indices.p)-1))-(sum(beta.indices.p)==1)),
                      1-as.numeric(beta.indices.s), rep(1,sum(1:(sum(beta.indices.s)-1))-(sum(beta.indices.s)==1)))
    
    # We now ensure that the previously selected main effects are kept and we include only potential pairwise interactions associated with previous main effects, but we allow for more main effects being selected still.
    pen.index <- c(rep(0,length(socios)), pen.index.ps)
    
    lasso_model_select <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, 
                                 penalty.factor = pen.index)
    
    AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
    minLambda <- lasso_model_select$lambda[which.min(AIC)]
    beta_select <- lasso_model_select$beta[,which.min(AIC)]
    
  }
  
  ## For extracting wald tests: We take the coxph output using the selected beta
  lp_exp <- paste(names(beta_select)[beta_select!=0],collapse = "+")
  if (lp_exp == ""){
    lp_exp <- 1
  }
  Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")",paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
  Routput <- coxph(Rformula,data=jsurvboot, x = FALSE, id = 1:nrow(jsurvboot), y = FALSE)
  
  Foutput <- list("coef" = coef(Routput), "vcov" = vcov(Routput))
  
  return(list(Foutput))
  
}

nonAssociationBootstrapShadow <- function(jsurvdata, socios, pollutants, event = "ASTHMA_H", P = 200){
  
  datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
  name.select <- pairwise.name.generator(pollutants)
  
  sampled_ids <- sample(1:nrow(datamatrix),replace=TRUE)
  
  boot_data <- datamatrix[sampled_ids,]
  shadowIds <- sample(1:nrow(boot_data),replace=FALSE)
  
  shadowPollutants <- boot_data[shadowIds,pollutants]
  shadows <- str_c("Shadow", 1:length(pollutants))
  colnames(shadowPollutants) <- shadows
  
  shadow.select <- pairwise.name.generator(shadows)
  
  boot_data <- cbind(boot_data, shadowPollutants)
  
  boot_data[,c(pollutants,shadows)] <- t((t(boot_data[,c(pollutants,shadows)] ) - colMeans(boot_data[,c(pollutants,shadows)]))/
                                           colSds(boot_data[,c(pollutants,shadows)] )) 
  
  jsurvboot <- jsurvdata[sampled_ids,]
  jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants
  jsurvboot <- cbind(jsurvboot, shadowPollutants) # And shadow pollutants
  
  strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))
  if (event == "ASTHMA_H"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)
  }
  
  if (event == "death"){
    jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$death)),strata = strata.indicator)
  }
  
  
  # New stuff with AIC instead:
  lasso_model_select <- glmnet(x=boot_data,y=jsurvStrat,family="cox",standardize=FALSE, 
                               penalty.factor = c(rep(0,length(socios)),rep(1,length(c(pollutants,shadows)))))
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  
  minLambdaMain <- lasso_model_select$lambda[which.min(AIC)]
  
  beta_select <- lasso_model_select$beta[,which.min(AIC)]
  beta_select_main_p <- beta_select[names(beta_select) %in% c(pollutants)]
  beta_select_main_s <- beta_select[names(beta_select) %in% c(shadows)]
  beta.indices.p <- (beta_select_main_p != 0)
  beta.indices.s <- (beta_select_main_s != 0)
  
  # SH
  
  Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants), "+",coalesce(pairwise.paster(pollutants[beta.indices.p]),"1"),"+", main.paster(shadows),"+",coalesce(pairwise.paster(shadows[beta.indices.s]),"1")),sep=""))
  MM_pol <- model.matrix(Pformula,jsurvboot)[,-1]
  MM_pol <- MM_pol[,order(match(colnames(MM_pol),c(name.select,shadow.select)))]
  MM <- cbind(boot_data[,1:length(socios)],MM_pol)
  
  pen.index.ps <- c(1-as.numeric(beta.indices.p), rep(1,sum(1:(sum(beta.indices.p)-1))-(sum(beta.indices.p)==1)),
                    1-as.numeric(beta.indices.s), rep(1,sum(1:(sum(beta.indices.s)-1))-(sum(beta.indices.s)==1)))
  
  # We now ensure that the previously selected main effects are kept and we include only potential pairwise interactions associated with previous main effects, but we allow for more main effects being selected still.
  pen.index <- c(rep(0,length(socios)), pen.index.ps)
  
  lasso_model_select <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, 
                               penalty.factor = pen.index)
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  minLambdaInt <- lasso_model_select$lambda[which.min(AIC)]
  
  # Bootstrap ranked lists when no association
  
  beta_select_NA <- matrix(nrow=P, ncol = length(name.select))
  
  Rout <- list()
  p <- 1
  while (p <= P){
    
    NA_lasso <- try(single.Perm.SH.Shadow(boot_data = boot_data, jsurvStrat = jsurvStrat, minLambdaMain = minLambdaMain, minLambdaInt = minLambdaInt, pollutants = pollutants, shadows = shadows, socios = socios, jsurvboot = jsurvboot, event = event), silent = TRUE)
    
    if(inherits(NA_lasso, "try-error")){
      next 
    }
    
    Rout[[p]] <- NA_lasso
    
    p <- p + 1
    
  }
  
  return(Rout) 
  
  
}


single.Perm.SH.Shadow <- function(boot_data,jsurvStrat, minLambdaMain, minLambdaInt, pollutants, shadows, socios, jsurvboot, event){
  
  name.select <- pairwise.name.generator(pollutants)
  shadow.select <- pairwise.name.generator(shadows)
  
  sample.indices <- sample(1:nrow(jsurvStrat),replace=FALSE)
  shadow.indices <- sample(1:nrow(jsurvStrat),replace=FALSE)
  
  NAstrat <- jsurvStrat[sample.indices,]
  
  jsurvbootNA <- jsurvboot
  jsurvbootNA$age_end_H <- jsurvboot$age_end_H[sample.indices]
  jsurvbootNA$SEX <- jsurvboot$SEX[sample.indices]
  jsurvbootNA$C_KOM <- jsurvboot$C_KOM[sample.indices]
  jsurvbootNA[,event] <- jsurvboot[sample.indices,event]
  
  jsurvbootNA[,shadows] <- jsurvboot[shadow.indices, shadows] # Three-way break each time between outcome, pollutants, and shadows
  
  #
  
  PformulaMain <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(c(pollutants,shadows))),sep=""))
  MMmain <- cbind(boot_data[,1:length(socios)],model.matrix(PformulaMain,jsurvbootNA)[,-1]) # Now important to use NA, since shadows are shuffled
  pen.index <- c(rep(0,length(socios)), rep(1,length(c(pollutants,shadows))))
  NA_lasso_Main <- glmnet(x=MMmain,y=NAstrat,family="cox",standardize=FALSE, 
                          penalty.factor = pen.index, lambda = minLambdaMain)
  
  beta_select <- NA_lasso_Main$beta[,1]
  beta_select_main_p <- beta_select[names(beta_select) %in% c(pollutants)]
  beta_select_main_s <- beta_select[names(beta_select) %in% c(shadows)]
  beta.indices.p <- (beta_select_main_p != 0)
  beta.indices.s <- (beta_select_main_s != 0)
  
  if (sum(c(beta.indices.p,beta.indices.s))>0){
    
    PformulaInt <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(c(pollutants,shadows)),"+",coalesce(pairwise.paster(pollutants[beta.indices.p]),"1"),"+",coalesce(pairwise.paster(shadows[beta.indices.s]),"1")),sep=""))
    
    MM_int_pol <- model.matrix(PformulaInt,jsurvbootNA)[,-1] # Now important to use jsurvbootNA, since shadows are shuffled
    MM_int_pol <- MM_int_pol[,order(match(colnames(MM_int_pol),c(name.select,shadow.select)))] # Reorder, so original pollutants go before shadows
    MMint <- cbind(boot_data[,1:length(socios)],MM_int_pol) # Socios are not shuffled, only outcome and shadows
    
    pen.index <- c(rep(0,length(socios)), c(1-as.numeric(beta.indices.p),rep(1,sum(1:(sum(beta.indices.p)-1))-(sum(beta.indices.p)==1))),
                   c(1-as.numeric(beta.indices.s),rep(1,sum(1:(sum(beta.indices.s)-1))-(sum(beta.indices.s)==1)))) 
    
    NA_lasso_int <- glmnet(x=MMint,y=NAstrat,family="cox",standardize=FALSE, 
                           penalty.factor = pen.index, lambda = minLambdaInt)
    
    ## For extracting coxph object:
    lp_exp <- paste(rownames(coef(NA_lasso_int))[coef(NA_lasso_int)[,1]!=0],collapse = "+", sep = "")
    if (lp_exp == ""){
      lp_exp <- 1
    }
    Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")", paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
    Routput <- coxph(Rformula,data=jsurvbootNA, x = FALSE, id = 1:nrow(jsurvbootNA), y = FALSE)
    
  }
  
  if (sum(c(beta.indices.p,beta.indices.s))==0){
    
    lp_exp <- paste(rownames(coef(NA_lasso_Main))[coef(NA_lasso_Main)[,1]!=0],collapse = "+", sep = "")
    if (lp_exp == ""){
      lp_exp <- 1
    }
    Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")", paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
    Routput <- coxph(Rformula,data=jsurvbootNA, x = FALSE, id = 1:nrow(jsurvbootNA), y = FALSE)
    
    
  }
  
  ## Return object;
  Foutput <- list("coef" = coef(Routput), "vcov" = vcov(Routput))
  
  #return(NA_lasso_int)
  return(Foutput)
}


rank.assignment.Shadow <- function(pollutants, shadows, SHbeta){ # input: pollutant names and list of model outputs
  
  ordersSH <- numeric(0) 
  pair.list <- pairwise.name.generator(pollutants)
  shadow.list <- pairwise.name.generator(shadows)
  
  for (i in (1:length(SHbeta))){
    
    beta <- SHbeta[[i]]$coef #coef(SHbeta[[i]])
    beta <- beta[names(beta) %in% c(pair.list, shadow.list)]
    
    VBeta <- SHbeta[[i]]$vcov #vcov(SHbeta[[i]])
    VBeta <- VBeta[rownames(VBeta) %in% c(pair.list, shadow.list), colnames(VBeta) %in% c(pair.list, shadow.list)]

    Fam.pvals <- rep(NA,length(names(beta)))
    names(Fam.pvals) <- names(beta)
    
    selected.mains <- names(Fam.pvals)[names(Fam.pvals) %in% c(pollutants,shadows)]
    for (j in selected.mains){
      
      k <- which(names(Fam.pvals)==j)
      
      R <- rep(0,length(beta))
      R[k] <- 1
      Q = 1
      
      tSize <- (R%*%beta)^2 * solve(R%*%(VBeta)%*%R)
      Fam.pvals[k] <- pchisq(tSize,df=Q,lower.tail = FALSE)
      
    }
    
    pairs.in.list <- names(Fam.pvals)[!names(Fam.pvals) %in% c(pollutants,shadows)]
    
    for(j in pairs.in.list){
      
      k <- which(names(Fam.pvals)==j)
      
      loc <- gregexpr(":",names(Fam.pvals)[k])[[1]][1]
      c1 <- substr(names(Fam.pvals)[k],1,loc-1)
      c2 <- substr(names(Fam.pvals)[k],loc+1,nchar(names(Fam.pvals)[k]))
      
      famIndices <- c(which(names(Fam.pvals) %in% c(c1,c2)), k)
      
      # Describing the hypothesis
      R <- matrix(rep(0,length(beta)),ncol = length(beta), nrow = length(famIndices))
      for (r in 1:length(famIndices)){
        R[r,famIndices[r]] <- 1
      }
      # Degrees of freedom
      Q = nrow(R)
      
      # Test statistic and p value
      tSize <- t((R%*%beta)) %*% solve(R%*%VBeta%*%t(R)) %*% (R%*%beta)
      Fam.pvals[k] <- pchisq(tSize, df = Q, lower.tail = FALSE) # df = length(famIndices)
      
    }
    
    # ... Now to fill in ranks to pollutants, with missings where we have not ranked a pollutant...
    # Following the original idea, now interactions represent the 3 member families, and the main effects represent 1 member families.
    
    full.pval.vec <- rep(1,length(c(pair.list, shadow.list)))
    names(full.pval.vec) <- c(pollutants, shadows, pair.list[!pair.list %in% pollutants], shadow.list[!shadow.list %in% shadows])
    
    full.pval.vec[names(full.pval.vec) %in% names(Fam.pvals)][order(match(full.pval.vec[names(full.pval.vec) %in% names(Fam.pvals)],names(Fam.pvals)))] <- Fam.pvals
    
    ord <- rank(full.pval.vec)
    ord[!names(ord) %in% names(Fam.pvals)] <- NA # Either not in contention or not selected when in contention!
    ordersSH <- rbind(ordersSH,ord)
    
  }
  
  return(ordersSH)
  
}
