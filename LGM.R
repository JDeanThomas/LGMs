require(OpenMx)
require(MASS)

# Means, Standard Deviations and Correlations from Duncan & Duncan (1996)

means <-  matrix(c(2.271, 2.560, 2.694, 2.965, .573, .554),1,6)
sds1  <-  c(1.002, .960, .912, .920, .504, .498)
Cor   <- matrix(c(1.000,  .640,  .586,  .454,  .001, -.214,			
                   .640, 1.000,  .670,  .566,  .038, -.149,
                   .586,  .670, 1.000,  .621,  .118, -.135,
                   .454,  .566,  .621, 1.000,  .091, -.163,
                   .001,  .038,  .118,  .091, 1.000, -.025,
                  -.214, -.149, -.135, -.163, -.025, 1.000), 6,6)

COV <- round(diag(sds1) %*% Cor %*% t(diag(sds1)),6)

data <- as.data.frame(mvrnorm(n = 321, mu = means, Sigma = COV, empirical = T))
colnames(data) <- c("year1", "year2", "year3", "year4", "Gender", "FamilyStatus" )
cov(data)



########
nVariables <- 4
nGrow      <- 2
vars  <-  c("year1", "year2", "year3", "year4")

# Specify psi - covariance of the latent growth parameters
PSIchol <- mxMatrix(type = "Lower", nrow = nGrow, ncol = nGrow, free = T, values = c(1,0,1), name = "PSIchol")
PSI     <- mxAlgebra(PSIchol%*% t(PSIchol), name = "PSI")
GrowCor <- mxAlgebra(cov2cor(PSI), name = "GrowCor")
# Specify the factor Loadings

growLoad <- c(1,1,1,1,
	          0,1,2,3)

LAMBDA  <- mxMatrix( type="Full", nrow=nVariables, ncol=nGrow, free=F, values = growLoad, name = "LAMBDA")  

# Specify the residuals of the manifest variables
EPSILON <- mxMatrix(type = "Diag", nrow = nVariables, ncol = nVariables, free = T, values = 1, name = "EPSILON")
expCov <- mxAlgebra(LAMBDA %&% PSI + EPSILON, name = "expCov")

# Specify the means of the latent factor and the manifest variables
GroMean <- mxMatrix(type="Full", nrow = nGrow, ncol = 1, free = T, values = c(1,.5), name = "GroMean")
ManMean  <- mxAlgebra(t(LAMBDA %*% GroMean), name = "ManMean")

# Data & Objective
dat       <- mxData(observed = data, type = "raw")                                                   # Data
GrowObj   <- mxFIMLObjective(covariance = "expCov", means ="ManMean", dimnames = vars)
             
growModel <- mxModel("LGM", PSIchol, PSI, GrowCor, LAMBDA, EPSILON, GroMean, ManMean, expCov, dat, GrowObj)
growFit  <- mxRun(growModel, unsafe = T)
summary(growFit)

mxEval(PSI, growFit)
mxEval(GrowCor, growFit)

#######################

nVariables <- 4
nGrow      <- 3
vars  <-  c("year1", "year2", "year3", "year4")

# Specify psi - covariance of the latent growth parameters
PSIchol <- mxMatrix(type = "Lower", nrow = nGrow, ncol = nGrow, free = T, values = c(1,0,0,1,0,1), name = "PSIchol")
PSI     <- mxAlgebra(PSIchol%*% t(PSIchol), name = "PSI")
GrowCor <- mxAlgebra(cov2cor(PSI), name = "GrowCor")
# Specify the factor Loadings

growLoad2 <- c(1,1,1,1,
	          -3,-1,1,3,
		       1,-1,-1,1)

LAMBDA  <- mxMatrix( type="Full", nrow=nVariables, ncol=nGrow, free=F, values = growLoad2, name = "LAMBDA")  

# Specify the residuals of the manifest variables
EPSILON <- mxMatrix(type = "Diag", nrow = nVariables, ncol = nVariables, free = T, values = 1, name = "EPSILON")
expCov <- mxAlgebra(LAMBDA %&% PSI + EPSILON, name = "expCov")

# Specify the means of the latent factor and the manifest variables
GroMean <- mxMatrix(type="Full", nrow = nGrow, ncol = 1, free = T, values = c(1,.5, .1), name = "GroMean")
ManMean  <- mxAlgebra(t(LAMBDA %*% GroMean), name = "ManMean")

# Data & Objective
dat       <- mxData(observed = data, type = "raw")                                                   # Data
GrowObj   <- mxFIMLObjective(covariance = "expCov", means ="ManMean", dimnames = vars)
             
growQModel <- mxModel("LGM", PSIchol, PSI, GrowCor, LAMBDA, EPSILON, GroMean, ManMean, expCov, dat, GrowObj)
growQFit  <- mxRun(growQModel, unsafe = T)
summary(growQFit)


mxEval(PSI, growQFit)
mxEval(GrowCor, growQFit)


# Plot the Results

#Grab the Means for the Latent Growth Parameters
means <- t(growQFit@output$matrices$LGM.GroMean)
# Grab the raw data
yearraw <- cbind(data$year1,data$year2,data$year3,data$year4)
# Create an expected mean vector
alc <- means[1] + means[2]*growLoad2[5:8] + means[3]*growLoad2[9:12]

# Graph the Effects

plot(growLoad2[5:8],alc, ylab = "Alcohol Use", xlab = "Year in High School", ylim = c(0,6), type = "b", xaxt = "n" , main = "Reanalysis of Duncan & Duncan (1996)")
axis(side = 1, at = growLoad2[5:8], labels = c("Freshman","Sophmore","Junior","Senior"))
for(i in 1:321){
	lines(growLoad2[5:8],yearraw[i,], ylab = " ", xlab = " ", ylim = c(0,2), type = "l", xaxt = "n", col = i , lwd = .5)
}
lines(growLoad2[5:8],alc, ylab = " ", xlab = " ", ylim = c(0,1), type = "b", xaxt = "n" , lwd = 3)



# LGC Model with Path Specification

nVariables<-6
nFactors<-2


require(OpenMx)

observed <- c("year1", "year2", "year3", "year4") #, "Gender", "FamilyStatus")
latent   <- c("intercept","slope")

pathData <- mxData(data, type="raw")

manResid <- mxPath(from=c(observed[1:4]), arrows=2,free=TRUE,values = c(1, 1, 1, 1),labels=c("res1","res2","res3","res4") )  # residual variances
GroCov   <- mxPath(from=c("intercept","slope"),arrows=2,connect="unique.pairs",free=TRUE,values=c(1, 1, 1),labels=c("varI", "cov", "varS"))     # latent variances and covariance
Int      <- mxPath(from="intercept",to=observed[1:4],arrows=1,free=FALSE,values=1 )                                                             # intercept loadings
Lin      <- mxPath(from="slope",to=observed[1:4],arrows=1,free=FALSE,values=c(0, 1, 2, 3))                                                   # slope loadings
ManMeans <- mxPath(from="one",to=observed[1:4] ,arrows=1,free=c(F,F,F,F),values=0 )                                                              # manifest means
LatMeans <- mxPath(from="one",to=c("intercept", "slope"),arrows=1,free=TRUE,values=c(1, 1),labels=c("meanI", "meanS") )                         # latent means


lgcModel <- mxModel("LGC",type="RAM", pathData, manifestVars= observed, latentVars= latent,
                    manResid, GroCov, Int, Lin, ManMeans, LatMeans) # model

lgcFit   <- mxRun(lgcModel)

summary(lgcFit)


### LGC Model with Predictors of the Growth Factors

observed <- c("year1", "year2", "year3", "year4", "Gender", "FamilyStatus")
latent   <- c("intercept","slope")

pathData <- mxData(data, type="raw")

manResid <- mxPath(from=c(observed), arrows=2,free=c(T,T,T,T,T,T),values = c(1, 1, 1, 1, 0, 0),labels=c("res1","res2","res3","res4","resG","resFS") )  # residual variances
GroCov   <- mxPath(from=c("intercept","slope"),arrows=2,connect="unique.pairs",free=TRUE,values=c(1, .1, 1),labels=c("varI", "cov", "varS"))     # latent variances and covariance
Int      <- mxPath(from="intercept",to=observed[1:4],arrows=1,free=FALSE,values=1 )                                                             # intercept loadings
Lin      <- mxPath(from="slope",to=observed[1:4],arrows=1,free=FALSE,values=c(0, 1, 2, 3))                                                   # slope loadings
ManMeans <- mxPath(from="one",to=observed ,arrows=1,free=c(F,F,F,F,T,T),values=0, labels = c("M1", "M2", "M3", "M4", "Ge", "FS") )                                                              # manifest means
LatMeans <- mxPath(from="one",to=latent,arrows=1,free=TRUE,values=c(1, .1),labels=c("meanI", "meanS") )                         # latent means

gender   <- mxPath(from="Gender", to = latent, arrows = 1, free = T, values = 0, labels = c("gen2int","gen2lin"))
family   <- mxPath(from="FamilyStatus", to = latent, arrows = 1, free = T, values = 0, labels = c("fam2int","fam2lin"))

lgc2Model <- mxModel("LGC",type="RAM", pathData, manifestVars= observed, latentVars= latent,
                    manResid, GroCov, Int, Lin, ManMeans, LatMeans, gender, family) # model

lgc2Fit   <- mxRun(lgc2Model)
summary(lgc2Fit)


### LGC Model with Predictors of the Growth Factors

observed <- c("year1", "year2", "year3", "year4")
latent   <- c("intercept","slope","quad")

pathData <- mxData(data, type="raw")

manResid <- mxPath(from=c(observed), arrows=2,free=c(T,T,T,T),values = c(1, 1, 1, 1),labels=c("res1","res2","res3","res4") )  # residual variances
GroCov   <- mxPath(from=latent,arrows=2,connect="unique.pairs",free=TRUE,values=.6 , labels = c("II", "IS","IQ","SS", "SQ","QQ"))     # latent variances and covariance
Int      <- mxPath(from="intercept",to=observed[1:4],arrows=1,free=FALSE,values=1 )                                                             # intercept loadings
Lin      <- mxPath(from="slope",to=observed[1:4],arrows=1,free=FALSE,values=c(-3, -1, 1, 3))                                                   # slope loadings
Quad     <- mxPath(from="quad",to=observed[1:4],arrows=1,free=FALSE,values=c(1, -1, -1, 1))                                                   # slope loadings
ManMeans <- mxPath(from="one",to=observed ,arrows=1,free=c(F,F,F,F),values=0, labels = c("M1", "M2", "M3", "M4") )                            # manifest means
LatMeans <- mxPath(from="one",to=latent,arrows=1,free=TRUE,values=c(1, .3,.1),labels=c("meanI", "meanS", "meanQ") )                         # latent means

lgc3Model <- mxModel("LGC",type="RAM", pathData, manifestVars= observed, latentVars= latent,
                    manResid, GroCov, Int, Lin, Quad, ManMeans, LatMeans) # model

lgc3Fit   <- mxRun(lgc3Model)
lgc3Fit   <- mxRun(lgc3Fit)

summary(lgc3Fit)

##$#$#$#$#$#$#$#$#
### LGC Model with Predictors of the Growth Factors

observed <- c("year1", "year2", "year3", "year4", "Gender", "FamilyStatus")
latent   <- c("intercept","slope","quad")

pathData <- mxData(data, type="raw")

manResid <- mxPath(from=c(observed), arrows=2,free=c(T,T,T,T,T,T),values = c(1, 1, 1, 1, 0, 0),labels=c("res1","res2","res3","res4","resG","resFS") )  # residual variances
GroCov   <- mxPath(from=latent,arrows=2,connect="unique.pairs",free=TRUE,values=.6 , labels = c("II", "IS","IQ","SS", "SQ","QQ"))     # latent variances and covariance
Int      <- mxPath(from="intercept",to=observed[1:4],arrows=1,free=FALSE,values=1 )                                                             # intercept loadings
Lin      <- mxPath(from="slope",to=observed[1:4],arrows=1,free=FALSE,values=c(-3, -1, 1, 3))                                                   # slope loadings
Quad     <- mxPath(from="quad",to=observed[1:4],arrows=1,free=FALSE,values=c(1, -1, -1, 1))                                                   # slope loadings


ManMeans <- mxPath(from="one",to=observed ,arrows=1,free=c(F,F,F,F,T,T),values=0, labels = c("M1", "M2", "M3", "M4", "Ge", "FS") )                                                              # manifest means
LatMeans <- mxPath(from="one",to=latent,arrows=1,free=TRUE,values=c(1, .3, .1),labels=c("meanI", "meanS", "meanQ") )                         # latent means

gender   <- mxPath(from="Gender", to = latent, arrows = 1, free = T, values = 0, labels = c("gen2int","gen2lin","gen2quad"))
family   <- mxPath(from="FamilyStatus", to = latent, arrows = 1, free = T, values = 0, labels = c("fam2int","fam2lin","fam2quad"))

lgcQModel <- mxModel("LGC",type="RAM", pathData, manifestVars= observed, latentVars= latent,
                    manResid, GroCov, Int, Lin, Quad, ManMeans, LatMeans, gender, family) # model


lgcQFit   <- mxRun(lgcQModel)
lgcQFit   <- mxRun(lgcQFit)

summary(lgcQFit)


##$#$#$#$#$#$#$#$

 # Latent Growth Model on 5 Repeatedly Measured Latent Traits
 nVariables <-20
 nFactors   <- 5
 nSubjects  <-500
 nGrow      <- 3
 
 # The covariance of y is lambda.y %*% solve(I-B)%*% (gamma%*%phi%*%gamma + psi) %*% t(solve(I-B))%*%t(lambda.y) + theta.epsilon
 
# Specify the lambda.y matrix - Factor loadings on the latent factors
lambda.y  <- matrix(c(1,rep(1,3), rep(0,16),                # You can change these values
                     rep(0,4) , 1,rep(1,3), rep(0,12),
                     rep(0,8) , 1,rep(1,3), rep(0,8) ,
                     rep(0,12), 1,rep(1,3), rep(0,4),
                     rep(0,16) , 1,rep(1,3))	,nrow=nVariables,ncol=nFactors)
 
# Specify gamma - Factor loadings on the exogenous variables (in this case growth parameters)
gamma  <- matrix(c(1,1,1,1,1,
 	              -2,-1,0,1,2,
 	               2,-1,-2,-1,2),5,3)

# Specify phi - covariance of the exogenous variables
phiChol <- matrix(c(1,.3,.15, 0,.5,.1,0,0,.1),3 )     # You can change these values
phi     <- phiChol%*% t(phiChol)

#Specify psi - residuals of the latent factors 
psi <- .4 * diag(5)                                    # You can change these values

# Specify theta.epsilon - the residuals of the manifest variables
theta.epsilon <- .5 * diag(nVariables)                # You can change these values

impCov <- lambda.y %*% (gamma%*%phi%*%t(gamma) + psi) %*% t(lambda.y) + theta.epsilon

 LatGrowMean <- matrix(c(2, -.5, .2),3,1)               # You can change these values
 LatFacMean <- gamma %*% LatGrowMean
 OBSMeans <- lambda.y %*% LatFacMean 
 
 
 data <- mvrnorm(1000, mu = OBSMeans, Sigma = impCov, empirical = T)
 colnames(data) <- paste(rep(paste("y", 1:5, sep = ""), each =4),1:4, sep='')
 


########
nVariables <- 20
nFactors   <- 5
nGrow      <- 3
selVars    <-  colnames(data)

loadF <- matrix(c(F,rep(T,3), rep(F,16),   rep(F,4) , F,rep(T,3), rep(F,12),   rep(F,8) , F,rep(T,3), rep(F,8) ,   
                    rep(F,12), F,rep(T,3), rep(F,4),   rep(F,16), F,rep(T,3)), nrow=nVariables,ncol=nFactors)
 
loadS <- matrix(c(1, rep(.3,3), rep(0,16),   rep(0,4) , 1, rep(.3,3), rep(0,12),   rep(0,8) , 1, rep(.3,3), rep(0,8) ,   
                     rep(0,12), 1, rep(.3,3), rep(0,4),   rep(0,16), 1, rep(.3,3)), nrow=nVariables, ncol=nFactors)

LAMBDA    <- mxMatrix(type="Full", nrow = nVariables, ncol = nFactors, free = loadF, values = loadS, name = "LAMBDA", dimnames = list(selVars, c(paste("F",1:5,sep = ""))))

# Specify gamma - Factor loadings on the exogenous variables (in this case growth parameters)
GAMMA    <- mxMatrix( type="Full", nrow=nFactors, ncol=nGrow, free=F, values=c(1,1,1,1,1,-2,-1,0,1,2,2,-1,-2,-1,2), name="GAMMA" ) # Matrix of Factor Loadings for the Growth Parameters


# Specify phi - covariance of the exogenous variables
PHIchol <- mxMatrix(type = "Lower", nrow = nGrow, ncol = nGrow, free = T, values = c(1,0,0,1,0,1), name = "PHIchol")
PHI     <- mxAlgebra(PHIchol%*% t(PHIchol), name = "PHI")

#Specify psi - residuals of the latent factors 
PSI <- mxMatrix(type = "Diag", nrow = nFactors, ncol = nFactors, free = T, values = 1, name = "PSI")



# Specify theta.epsilon - the residuals of the manifest variables
EPSILON <- mxMatrix(type = "Diag", nrow = nVariables, ncol = nVariables, free = T, values = 1, name = "EPSILON")


expCov <- mxAlgebra(LAMBDA %*% (GAMMA %*% PHI %*% t(GAMMA) + PSI) %*%  t(LAMBDA) + EPSILON, name = "expCov")


GroMean <- mxMatrix(type="Full", nrow = nGrow, ncol = 1, free = T, values = c(1,.5,.1), name = "GroMean")
FacMean  <- mxAlgebra(GAMMA %*% GroMean, name = "FacMean", dimnames = list(paste("F",1:5, sep=""),"FacMeans"))
ManMean  <- mxAlgebra(t(LAMBDA %*% FacMean), "ManMean")

           # Data & Objective
dat       <- mxData(observed = data, type = "raw")                                                                               # Data
GrowObj   <- mxFIMLObjective(covariance = "expCov", means ="ManMean", dimnames = colnames(data))                                 # Objective 


facgrowModel <- mxModel("LGM", LAMBDA, GAMMA, PHIchol, PHI, PSI, EPSILON,
                            ManMean, GroMean, FacMean, expCov, dat, GrowObj)

facgrowFit  <- mxRun(facgrowModel, unsafe = T)
#facgrowFit  <- mxRun(facgrowFit, unsafe = T)

summary(facgrowFit)



##$#$#$#$#$#$#$#$
 # Latent Growth Simplex Model on 5 Repeatedly Measured Latent Traits
 nVariables <-20
 nFactors   <- 5
 nSubjects  <-500
 nGrow      <- 3
 
 # The covariance of y is lambda.y %*% solve(I-B)%*% (gamma%*%phi%*%gamma + psi) %*% t(solve(I-B))%*%t(lambda.y) + theta.epsilon
 
# Specify the lambda.y matrix - Factor loadings on the latent factors
lambda.y  <- matrix(c(1,rep(1,3), rep(0,16),                # You can change these values
                     rep(0,4) , 1,rep(1,3), rep(0,12),
                     rep(0,8) , 1,rep(1,3), rep(0,8) ,
                     rep(0,12), 1,rep(1,3), rep(0,4),
                     rep(0,16) , 1,rep(1,3))	,nrow=nVariables,ncol=nFactors)
 
# Specify the I-B^-1  - The causal pathways between the latent factors

I5 <- diag(5)
B <- matrix(c(0,.5,0,0,0,                             # You can change these values
 	          0,0,.5,0,0,
 	          0,0,0,.5,0,
 	          0,0,0,0,.5,
 		      0,0,0,0,0),5)
 
# Specify gamma - Factor loadings on the exogenous variables (in this case growth parameters)
gamma  <- matrix(c(1,1,1,1,1,
 	              -2,-1,0,1,2,
 	               2,-1,-2,-1,2),5,3)

# Specify phi - covariance of the exogenous variables
phiChol <- matrix(c(1,.3,.15, 0,.5,.1,0,0,.1),3 )     # You can change these values
phi     <- phiChol%*% t(phiChol)

#Specify psi - residuals of the latent factors 
psi <- .4 * diag(5)                                    # You can change these values

# Specify theta.epsilon - the residuals of the manifest variables
theta.epsilon <- .5 * diag(nVariables)                # You can change these values

impCov <- lambda.y %*% solve(I5-B)%*% (gamma%*%phi%*%t(gamma) + psi) %*% t(solve(I5-B))%*%t(lambda.y) + theta.epsilon

 LatGrowMean <- matrix(c(2, -.5, .2),3,1)               # You can change these values
 LatFacMean <- solve(I5-B)%*%gamma %*% LatGrowMean
 OBSMeans <- lambda.y %*% LatFacMean 
 
 
 data <- mvrnorm(1000, mu = OBSMeans, Sigma = impCov, empirical = T)
 colnames(data) <- paste(rep(paste("y", 1:5, sep = ""), each =4),1:4, sep='')
 


########



nVariables <- 20
nFactors   <- 5
nGrow      <- 3
selVars    <-  colnames(data)

loadF <- matrix(c(F,rep(T,3), rep(F,16),   rep(F,4) , F,rep(T,3), rep(F,12),   rep(F,8) , F,rep(T,3), rep(F,8) ,   
                    rep(F,12), F,rep(T,3), rep(F,4),   rep(F,16), F,rep(T,3)), nrow=nVariables,ncol=nFactors)
 
loadS <- matrix(c(1, rep(.3,3), rep(0,16),   rep(0,4) , 1, rep(.3,3), rep(0,12),   rep(0,8) , 1, rep(.3,3), rep(0,8) ,   
                     rep(0,12), 1, rep(.3,3), rep(0,4),   rep(0,16), 1, rep(.3,3)), nrow=nVariables, ncol=nFactors)

LAMBDA    <- mxMatrix(type="Full", nrow = nVariables, ncol = nFactors, free = loadF, values = loadS, name = "LAMBDA", dimnames = list(selVars, c(paste("F",1:5,sep = ""))))
IDEN      <- mxMatrix(type="Iden", nrow = nFactors, ncol = nFactors, name = "IDEN")

betaF   <- matrix(c(F,T,F,F,F,   F,F,T,F,F,   F,F,F,T,F,   F,F,F,F,T,   F,F,F,F,F),5)
betaS   <- matrix(c(0,.1,0,0,0,   0,0,.1,0,0,   0,0,0,.1,0,   0,0,0,0,.1,   0,0,0,0,0),5)
betalab <- matrix(c(NA,"AR1",NA,NA,NA,   NA,NA,"AR2",NA,NA,   NA,NA,NA,"AR3",NA,   NA,NA,NA,NA,"AR4",   NA,NA,NA,NA,NA),5)

								   
BETA     <- mxMatrix(type="Full", nrow = nFactors, ncol = nFactors, free = betaF, values = betaS, name = "BETA", ubound = 1, labels = betalab)

# Specify gamma - Factor loadings on the exogenous variables (in this case growth parameters)
GAMMA    <- mxMatrix( type="Full", nrow=nFactors, ncol=nGrow, free=F, values=c(1,1,1,1,1,-2,-1,0,1,2,2,-1,-2,-1,2), name="GAMMA" ) # Matrix of Factor Loadings for the Growth Parameters


# Specify phi - covariance of the exogenous variables
PHIchol <- mxMatrix(type = "Lower", nrow = nGrow, ncol = nGrow, free = T, values = c(1,0,0,1,0,1), name = "PHIchol")
PHI     <- mxAlgebra(PHIchol%*% t(PHIchol), name = "PHI")

#Specify psi - residuals of the latent factors 
PSI <- mxMatrix(type = "Diag", nrow = nFactors, ncol = nFactors, free = T, values = 1, name = "PSI")



# Specify theta.epsilon - the residuals of the manifest variables
EPSILON <- mxMatrix(type = "Diag", nrow = nVariables, ncol = nVariables, free = T, values = 1, name = "EPSILON")


expCov <- mxAlgebra(LAMBDA %*% solve(IDEN-BETA)%*% (GAMMA %*% PHI %*% t(GAMMA) + PSI) %*% t(solve(IDEN-BETA))%*%t(LAMBDA) + EPSILON, name = "expCov")


GroMean <- mxMatrix(type="Full", nrow = nGrow, ncol = 1, free = T, values = c(1,.5,.1), name = "GroMean")
FacMean  <- mxAlgebra(solve(IDEN-BETA) %*% GAMMA %*% GroMean, name = "FacMean", dimnames = list(paste("F",1:5, sep=""),"FacMeans"))
ManMean  <- mxAlgebra(t(LAMBDA %*% FacMean), "ManMean")

           # Data & Objective
dat       <- mxData(observed = data, type = "raw")                                                                               # Data
GrowObj   <- mxFIMLObjective(covariance = "expCov", means ="ManMean", dimnames = colnames(data))                                 # Objective 


growModel <- mxModel("LGM", LAMBDA, BETA, GAMMA, PHIchol, PHI, PSI, EPSILON, IDEN,
                            ManMean, GroMean, FacMean, expCov, dat, GrowObj)

growFit  <- mxRun(growModel, unsafe = T)
#growFit  <- mxRun(growFit, unsafe = T)

summary(growFit)


