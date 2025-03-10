# This script contains dependencies and functions that are necessary for doing entropy based rank agreement selection.
# Towards the end of the script we compute a selection set with the entropy rank agreement selection method.
# For illustration we do the same with the sequential rank agreement and we give examples of how to use the free lasso and the hierarchical lasso on the same data.

library(survival)
library(Matrix)
library(matrixStats)
library(glmnet)
library(dplyr)
library(stringr)
library(rlang)
library(entropy)
library(SuperRanker)
library(mvtnorm)

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

residual.rank.randomizer <- function(orders){
  
  for (i in (1:nrow(orders))[rowSums(is.na(orders))==1]){
    
    orders[i,is.na(orders[i,])] <- ncol(orders)
    
  }
  
  for (i in (1:nrow(orders))[rowSums(is.na(orders))>0]){
    
    orders[i,is.na(orders[i,])] <- sample((sum(!is.na(orders[i,]))+1):ncol(orders),replace=FALSE)
    
  }
  
  return(orders)
  
}

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
  form <- ""
  
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
  form <- ""
  
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
  
  if (sum(beta.indices)>0){
    
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
  
  if (sum(beta.indices)>0){
    
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
  
  if (sum(beta.indices)==0){
    
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


# Artifical data example:

set.seed(1516792)

# Setup:

N <- 5000
nTry <- 200
p <- 10
r <- 0.95
M <- 200
sigma <- matrix(nrow=p,ncol=p)

for (i in 1:nrow(sigma)){
  
  sigma[i,i] <- 1
  sigma[i,-i] <- r
  
}

X <- rmvnorm(n=N,mean = rep(2,nrow(sigma)), sigma = sigma)
X[which(X<0, arr.ind=TRUE)] <- 0

SEX <- rbinom(n=N,size=1,prob=0.5)
h <- 0.1 + 0.1*X[,1]+0.2*X[,2]+0.1*X[,1]*X[,2] + 0.1*SEX
S <- rexp(N,h)
H <- 0.6
Y <- as.numeric(S<=H)
C <- (1 - Y)

jsurvdata = data.frame("ASTHMA_H" = Y, "age_end_H" = rowMins(cbind(S,0.6)), SEX = SEX, C_KOM = 1, X)
socios <- c("SEX")
pollutants <- str_c("X",1:10)
event <- "ASTHMA_H"


# The entropy rank agreement method:

# We carry out invididual selections:

try.result <- iLassoBoot(jsurvdata, socios = socios, pollutants = pollutants, event = "ASTHMA_H")

for (i in 2:nTry){
  
  try.result.add <- iLassoBoot(jsurvdata, socios = socios, pollutants = pollutants, event = "ASTHMA_H")
  
  try.result <- append(try.result, try.result.add)
  
}

try.NA <- list()

for (i in 1:nTry){
  
  try.NA[[i]] <- nonAssociationBootstrap(jsurvdata, socios = socios, pollutants = pollutants, event = "ASTHMA_H")
  
}

# We compute orders lists from the generated selections:

orders.result <- rank.assignment(pollutants, try.result)

orders.NA <- list()

for (i in 1:length(try.NA)){
  
  orders.NA[[i]] <- rank.assignment(pollutants, try.NA[[i]])
  
}

# Now for the depth estimation part:

e.result <- rankEntropy(d = ncol(orders.result), ranks = orders.result, P = 1000)

e.NA <- matrix(nrow = length(orders.NA), ncol = ncol(orders.result))

for (i in 1:nrow(e.NA)){
  
  e.NA[i,] <- rankEntropy(d = ncol(orders.NA[[i]]), ranks = orders.NA[[i]], P = 1000)
  
}

e.NA.mean <- colMeans(e.NA)
e.NA.sd <- colSds(e.NA)
e.NA.dl <- e.NA.mean - 1.96*e.NA.sd

d.raw.era <- which(e.result > e.NA.dl)[1] - 1


# Based on the depth estimate we select the final predictors:

if (length(d.raw.era) == 0) { d.raw.era = length(e.result) }
if (is.na(d.raw.era)) { d.raw.era = length(e.result)}

d.era <- d.raw.era

if (d.raw.era >0){
  
  d <- d.raw.era
  
  indices <- matrix(nrow = nrow(orders.result), ncol = d)
  coef.names = pairwise.name.generator(pollutants)
  coef.selection <- numeric(0)
  selectionFrame <- numeric(0)
  
  for (i in 1:nrow(orders.result)){
    
    if (sum(!is.na(orders.result[i,]))>0){
      
      addRow <-  coef.names[which(orders.result[i,]<=d)]
      coef.selection <- c(coef.selection, addRow)
      
      if (length(addRow)<d){
        
        addRow <- c(addRow, rep(NA, d-length(addRow)))
        
      }
      
      selectionFrame <- rbind(selectionFrame, addRow)
    }
    
  }
  
  table(coef.selection)
  threshold.pollutants <- table(coef.selection)[order(table(coef.selection),decreasing=T)][1:d]
  marginal.era <- threshold.pollutants
  
  # Translation from pollutants to families in threshold methods:
  
  hierarchical.selection <- marginal.era 
  
  addOns <- numeric(0)
  
  for (h in 1:length(hierarchical.selection)){
    
    addOns <- c(addOns, pollutants[which(string.contains(pollutants,hierarchical.selection[h])>0)])
    
  }
  
  addOns <- unique(addOns)
  addOns <- addOns[!addOns %in% hierarchical.selection]
  
  hierarchical.selection <- c(hierarchical.selection, addOns)
  
  marginal.era <- names(hierarchical.selection)
  
}

# The object marginal.era now contains the selected predictors with entropy based variable selection.


# Now with the sequential rank agreement instead of the entropy rank agreement:
# This is a different way of estimating the depth.

sra.result <- sra(object = t(orders.result), B = 1000)

sra.NA <- matrix(nrow = length(orders.NA), ncol = ncol(orders.result))

for (i in 1:nrow(sra.NA)){
  
  sra.NA[i,] <- sra(object = t(orders.NA[[i]]), B = 1000) #rankEntropy(d = ncol(orders.NA[[i]]), ranks = orders.NA[[i]], P = 1000)
  
}

sra.NA.mean <- colMeans(sra.NA)
sra.NA.sd <- colSds(sra.NA)
sra.NA.dl <- sra.NA.mean - 1.96*sra.NA.sd

d.raw.sra <- which(sra.result > sra.NA.dl)[1] - 1

# With the sequential rank agreement depth estimate we construct the predictor selection set:

if (length(d.raw.sra) == 0){ d.raw.sra = length(sra.result)}
if (is.na(d.raw.sra)) { d.raw.sra = length(sra.result)}

d.sra <- d.raw.sra

if (d.raw.sra > 0){
  
  d <- d.raw.sra
  
  indices <- matrix(nrow = nrow(orders.result), ncol = d)
  coef.names = pairwise.name.generator(pollutants)
  coef.selection <- numeric(0)
  selectionFrame <- numeric(0)
  
  for (i in 1:nrow(orders.result)){
    
    if (sum(!is.na(orders.result[i,]))>0){
      
      addRow <-  coef.names[which(orders.result[i,]<=d)]
      coef.selection <- c(coef.selection, addRow)
      
      if (length(addRow)<d){
        
        addRow <- c(addRow, rep(NA, d-length(addRow)))
        
      }
      
      selectionFrame <- rbind(selectionFrame, addRow)
    }
    
  }
  
  threshold.pollutants <- table(coef.selection)[order(table(coef.selection),decreasing=T)][1:d]
  marginal.sra <- threshold.pollutants
  
  # Translation from pollutants to families in threshold methods:
  
  hierarchical.selection <- marginal.sra 
  
  addOns <- numeric(0)
  
  for (h in 1:length(hierarchical.selection)){
    
    addOns <- c(addOns, pollutants[which(string.contains(pollutants,hierarchical.selection[h])>0)])
    
  }
  
  addOns <- unique(addOns)
  addOns <- addOns[!addOns %in% hierarchical.selection]
  
  hierarchical.selection <- c(hierarchical.selection, addOns)
  
  marginal.sra <- names(hierarchical.selection)
  
}

# The object marginal.sra now contains the selected predictors with variance based variable selection.


# Application of basic methods:

# A single free lasso:

datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
name.select <- pairwise.name.generator(pollutants)
pair.list <- pairwise.paster(pollutants)

boot_data <- datamatrix
jsurvboot <- jsurvdata
boot_data[,pollutants] <- t((t(boot_data[,pollutants] ) - colMeans(boot_data[,pollutants] ))/colSds(boot_data[,pollutants] )) # standardizing a priori, as to not standardize the interaction terms themselves (thereby introducing varying scales), but only the separate terms.
jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants

strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))

jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)

Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants),"+",coalesce(pairwise.paster(pollutants),"1")),sep=""))
MM <- cbind(boot_data[,1:length(socios)],model.matrix(Pformula,jsurvboot)[,-1])

pen.index <- c(rep(0,length(socios)),rep(1,length(c(pairwise.name.generator(pollutants))))) 

lasso_model <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = pen.index)

AIC <- deviance(lasso_model) + 2 * lasso_model$df

beta_select_free <- lasso_model$beta[,which.min(AIC)]
beta_select_free <- names(beta_select_free)[names(beta_select_free) %in% pairwise.name.generator(pollutants)]


# A single hierarchical lasso:

datamatrix <- as.matrix(jsurvdata[,c(socios,pollutants)])
name.select <- pairwise.name.generator(pollutants)
pair.list <- pairwise.paster(pollutants)

boot_data <- datamatrix
boot_data[,pollutants] <- t((t(boot_data[,pollutants] ) - colMeans(boot_data[,pollutants] ))/colSds(boot_data[,pollutants] )) # standardizing a priori, as to not standardize the interaction terms themselves (thereby introducing varying scales), but only the separate terms.
jsurvboot <- jsurvdata
jsurvboot[,pollutants] <- boot_data[,pollutants] # Taking the standardized pollutants

strata.indicator <- as.numeric(strata(jsurvboot$SEX,jsurvboot$C_KOM))

jsurvStrat <- stratifySurv(Surv(as.numeric(jsurvboot$age_end_H),as.numeric(jsurvboot$ASTHMA_H)),strata = strata.indicator)

lasso_model_select <- glmnet(x=boot_data,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = c(rep(0,length(socios)),rep(1,length(pollutants))))

AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df

minLambda <- lasso_model_select$lambda[which.min(AIC)]
beta_select <- lasso_model_select$beta[,which.min(AIC)]
beta_select_main <- beta_select[names(beta_select) %in% pollutants]
beta.indices <- (beta_select_main != 0)

if (sum(beta.indices)>0){
  
  Pformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")~",main.paster(pollutants),"+",coalesce(pairwise.paster(pollutants[beta.indices]),"1")),sep=""))
  MM <- cbind(boot_data[,1:length(socios)],model.matrix(Pformula,jsurvboot)[,-1])
  
  pen.index <- c(rep(0,length(socios)),c(1-as.numeric(beta.indices),rep(1,sum(1:(sum(beta.indices)-1))-(sum(beta.indices)==1)))) # Removing the 1 for no. interactions that happens when exactly one pollutant is chosen!
  
  lasso_model_select <- glmnet(x=MM,y=jsurvStrat,family="cox",standardize=FALSE, penalty.factor = pen.index)
  
  AIC <- deviance(lasso_model_select) + 2 * lasso_model_select$df
  minLambda <- lasso_model_select$lambda[which.min(AIC)]
  beta_select <- lasso_model_select$beta[,which.min(AIC)]
  
}

lp_exp <- paste(names(beta_select)[beta_select!=0],collapse = "+")
if (lp_exp == ""){
  lp_exp <- 1
}
Rformula <- as.formula(paste(str_c("Surv(age_end_H,",event,")",paste("~strata(SEX,C_KOM)",lp_exp, sep = "+")),sep=""))
Routput <- coxph(Rformula,data=jsurvboot, x = FALSE, id = 1:nrow(jsurvboot), y = FALSE)

Foutput <- list("coef" = coef(Routput), "vcov" = vcov(Routput))

beta_select_hierarchical <- names(Foutput$coef)[names(Foutput$coef) %in% pairwise.name.generator(pollutants)]

