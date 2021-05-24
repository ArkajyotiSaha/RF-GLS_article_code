rm(list = ls())
list.of.packages <- c("Rcpp", "RcppArmadillo", "data.table", "matrixStats", "BRISC", "lhs", "rdist", "MASS", "emdbook", "RandomForestsGLS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rcpp)
library(RcppArmadillo)
library(data.table)
setDTthreads(1)
Rcpp::sourceCpp("final_spRF_new_bigger_oral.cpp")

#temp <- commandArgs(TRUE)
#iter_tot <- (as.numeric(temp[1])-1)*24+as.numeric(temp[2])
iter_tot <- 1 #1..2700
indica <- iter_tot

fried <- function(xx)
{
  y <- rep(0, nrow(xx))
  for(i in 1:nrow(xx)){
    x1 <- xx[i,1]
    x2 <- xx[i,2]
    x3 <- xx[i,3]
    x4 <- xx[i,4]
    x5 <- xx[i,5]
    
    x6 <- xx[i,6]
    x7 <- xx[i,7]
    x8 <- xx[i,8]
    x9 <- xx[i,9]
    x10 <- xx[i,10]
    
    x11 <- xx[i,11]
    x12 <- xx[i,12]
    x13 <- xx[i,13]
    x14 <- xx[i,14]
    x15 <- xx[i,15]
    
    
    term1 <- 10 * sin(pi*x1*x2)
    term2 <- 20 * (x3-0.5)^2
    term3 <- 10*x4
    term4 <- 5*x5
    term5 <- 3 * 1/((x6+1)*(x7+1))
    term6 <- 4 * exp(x8^2)
    term7 <- 30 * x9^2 * x10
    
    term8 <- 5 * (exp(x11^2) * sin(pi*x12) + exp(x12^2) * sin(pi*x11))
    term9 <- 10 * x13^2 * cos(pi*x14)
    term10 <- 20 * x15^4
    
    
    y[i] <- (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10)/6
  }
  return(y)
}

Coeffmatrix<-expand.grid(iter = 1:100, sigma.sq = c(1,5,10), fraction = c(1/4,2/4,3/4),
                         percent = c(1/100,10/100,25/100))

iter <- Coeffmatrix$iter[indica]
n = 250


set.seed(iter+1)
library(lhs)
xmat <- randomLHS(n, 15)


set.seed(iter+2)
s1 <- runif(n,0, 1)
set.seed(iter+3)
s2 <- runif(n,0, 1)

smat <- cbind(s1, s2)

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}


sigma.sq = Coeffmatrix$sigma.sq[indica]
phi = 3/(sqrt(2)*Coeffmatrix$fraction[indica])
tau.sq = Coeffmatrix$sigma.sq[indica] * Coeffmatrix$percent[indica]
D200 <- as.matrix(dist(smat))
nu <- 3/2
R200 <- (D200*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D200*phi, nu=nu)
diag(R200) <- 1
w200 <- (sigma.sq)^0.5 * rmvn(1, rep(0,n), R200)

set.seed(iter+4)
y200 <- rnorm(n, fried(xmat) + w200, sqrt(tau.sq))


set.seed(iter+101)
simp_row <- sample(1:10)
set.seed(iter+102)
simp_col <- sample(1:10)

stest_index <- 0

for(ip in 1:10){
  stest_index <- c(stest_index, which(smat[,1]< simp_row[ip]*0.1 & smat[,1] >= (simp_row[ip]-1)*0.1 & smat[,2]< simp_col[ip]*0.1  & smat[,2]>= (simp_col[ip]-1)*0.1))
}

stest_index <- stest_index[-1]
stest <- smat[stest_index,]
strain <- smat[-stest_index,]

xtrain200 <- xmat[-stest_index,]
strain200 <- strain
ytrain200 <- y200[-stest_index]
fxtrain200  <- fried(xtrain200)

xtest200 <- xmat[stest_index,]
stest200 <- stest
ytest200 <- y200[stest_index]
fxtest200 <- fried(xtest200)

set.seed(iter+5)
xestimation200 <- randomLHS(1000, 15)
colnames(xestimation200) <- colnames(xtrain200)
fxestimation200 <- fried(xestimation200)



ntree <- 100
nthsize <- 5
n.neighbors.number <- 20

set.seed(iter+7)
sp <- randomForest::randomForest(xtrain200, ytrain200, ntree = ntree,nodesize = nthsize)

rftest200 <- predict(sp, xtest200)
rfestimation200 <- predict(sp, xestimation200)
rf_estimation200_res <- mean((rfestimation200 - fxestimation200)^2, trim = 0.0) #rf_estimation_error
rf_res <- mean((rftest200 - ytest200)^2, trim = 0.0)

rftrain200 <- predict(sp, xtrain200)
rf_residual <- ytrain200 - rftrain200

library(rdist)
Dis <- pdist(strain200)
dis_pred <- cdist(stest200, strain200)
set.seed(iter+7)

sp_hengl <- randomForest::randomForest(cbind(xtrain200,Dis), ytrain200, xtest = cbind(xtest200,dis_pred), ntree = ntree,nthsize=nthsize, importance = T, localImp=T)
#sp_out_hengl <- predict(sp_hengl, cbind(xtest200,dis_pred))
rf_res_hengl <- mean((sp_hengl$test$predicted - ytest200)^2, trim = 0.0)
sp_out_hengl <- sp_hengl$test$predicted

set.seed(iter+7)
sp_lat_long <- randomForest::randomForest(cbind(xtrain200,strain200), ytrain200, xtest = cbind(xtest200, stest200), ntree = ntree,nodesize = nthsize, importance = T, localImp= T)
rf_res_lat_long <- mean((sp_lat_long$test$predicted - ytest200)^2, trim = 0.0)
sp_out_lat_long <- sp_lat_long$test$predicted


library(BRISC)
est_theta <- BRISC_estimation(strain200, y = c(rf_residual), x = matrix(1, nrow(xtrain200), 1), n.neighbors = n.neighbors.number)

total_x_prediction <- rbind(xtest200, xestimation200)


library(matrixStats)
spatialrandomforest_tot <- function(x,y, coords, xtest, nsample, mdim, ntest, nthsize = 5, ntree = 50, mtry = 1,sigmasq = 1, tausq = 0.1, phi = 1000, nu = 0.5, covModel = 0, ytest, fxin, fxout, trim = 0.02){
  nrnodes <- 2*nsample + 1
  rfull <- regRFtestnew(x, y, mdim, nsample, nthsize, nrnodes, ntree, mtry, sigmasq, tausq, phi, nu, covModel, coords, xtest, ntest)
  list(prediction.in.mse.median = mean((rowMedians(rfull[[1]]) - fxin)^2, trim = trim), prediction.in.mse.mean = mean((rowMeans(rfull[[1]]) - fxin)^2, trim = trim), prediction.out.mse.medians = mean((rowMedians(rfull[[2]]) - fxout)^2, trim = trim), prediction.out.mse.mean = mean((rowMeans(rfull[[2]]) - fxout)^2, trim = trim), rfull[[2]], rfull[[1]])
}

set.seed(iter+8)
sprf_res <- spatialrandomforest_tot(x = xtrain200,y = ytrain200, coords = strain200, xtest = total_x_prediction, nsample = nrow(xtrain200), nthsize = nthsize,
                                    mdim = ncol(xtrain200), ntree = ntree, ntest = nrow(total_x_prediction), mtry = floor(ncol(xtrain200)/3),sigmasq = est_theta$Theta[1], tausq = est_theta$Theta[2], phi = est_theta$Theta[3], ytest = NULL, fxin = fxtrain200, fxout = c(fxtest200, fxestimation200), trim = 0.00)


rfgls_unknown_test200_res <- mean((rowMeans(sprf_res[[5]])[1:length(fxtest200)] - fxtest200)^2, trim = 0.0)
rfgls_unknown_estimation200_res <- mean((rowMeans(sprf_res[[5]])[(length(fxtest200)+1):nrow(sprf_res[[5]])] - fxestimation200)^2, trim = 0.0) #rfgls_estimation_error

est_theta_ols <- est_theta

rfgls_insample <- rowMeans(sprf_res[[6]])
rfgls_residual <- ytrain200 - rfgls_insample

est_theta_gls <- BRISC_estimation(strain200, y = c(rfgls_residual), x = matrix(1, nrow(xtrain200), 1), n.neighbors = n.neighbors.number)


pred_gls_spatial <- BRISC_prediction(BRISC_Out = est_theta_gls, X.0 = matrix(1, nrow(xtest200), 1), coords.0 = stest200)
pred_ols_spatial <- BRISC_prediction(BRISC_Out = est_theta_ols, X.0 = matrix(1, nrow(xtest200), 1), coords.0 = stest200)

pred_ols = rftest200 + pred_ols_spatial$prediction
pred_gls = rowMeans(sprf_res[[5]])[1:length(fxtest200)] + pred_gls_spatial$prediction

rf_res_spatial <- mean((pred_ols - ytest200)^2) #rf_rk_prediction_error
rfgls_res_spatial <- mean((pred_gls - ytest200)^2) #rfgls_prediction_error
