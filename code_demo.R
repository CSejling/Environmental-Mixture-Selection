# In this file we review how to run some pieces of code that are the building blocks in some of the analyses of the paper...

source("generative_code.R")

# ------- #

# Run of a competing risk cox model: The socios are the socioeconomic counfounding factors to adjust for and pollutants are the ambient air pollutants to include in the model

out <- CSC(list(as.formula(paste(str_c("Hist(follow_up_time,status)~strata(SEX,C_KOM)+",paste(socios,collapse="+"),"+",paste(pollutants, collapse ="+")),sep="")),
                as.formula(paste(str_c("Hist(follow_up_time,status)~strata(SEX,C_KOM)+",paste(socios,collapse="+"),"+",paste(pollutants, collapse ="+")),sep=""))),
           surv.type = "hazard", data = ...)

# To predict absolute risks:

predictRisk(object = out, newdata = ..., cause = 1, times = c(5,10,15))


# ------- #


# Run of single hierarchical lasso:

# jsurvdata is the data sample, socios is a list of variable names corresponding with confounders to be adjusted for, and pollutants is the list of pollutants to include with interactions in the model
try(iLassoBoot(jsurvdata=..., socios = ..., pollutants = ...,  lambda.list = NULL), silent = TRUE)

# ranking a sample of iLasso outputs: SHbeta is a list of coefficient estimates from iLasso runs.

rank.assignment(pollutants = ..., SHbeta = ...)


# --------- #


# Calculation of entropy rank agreement: the ranks is a matrix of ranks as outputted by rank.assignment, d is the depth to which we want to go, and P is the number of residual rank randomizations

rankEntropy(d=..., ranks = ..., P = ...)


# --------- #


# Shadow copy analogues: To deal with the shadow copy, we need to consider some small adjustments.

# The shadow version of the iLasso runs the same. Inside the function, a shadow copy is generated
try(iLassoBootShadow(jsurvdata=..., socios = ..., pollutants = ...), silent = TRUE)

# For the rank assignment we need to specify which predictors are shadows and name them correctly. This is to not confuse families of pollutants with each other, when there are shadows present, since shadows should never interact with real pollutants.
shadows <- str_c("Shadow", 1:length(pollutants))
ordersSH <- rank.assignment.Shadow(pollutants = pollutants, shadows = shadows, SHbeta = SHbeta)

