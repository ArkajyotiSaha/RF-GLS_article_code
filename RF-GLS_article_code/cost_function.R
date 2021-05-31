rfgls_cost_t1 <- function(c_0,y,x, rho, sigma.sq){
  n <- length(y)
  z <- as.numeric(x > c_0)
  zQz <- matrix(0,2,2)
  zQz[1,1] <- xQy(z, z, rho, sigma.sq)
  zQz[1,2] <- xQy(z, 1-z, rho, sigma.sq)
  zQz[2,1] <- xQy(1-z, z, rho, sigma.sq)
  zQz[2,2] <- xQy(1-z, 1-z, rho, sigma.sq)
  zQy <- matrix(0,2,1)
  zQy[1,1] <- xQy(z, y, rho, sigma.sq)
  zQy[2,1] <- xQy(1-z, y, rho, sigma.sq)
  beta <- solve(zQz, zQy)
  cost <- xQy(y - cbind(z,1-z) %*% beta, y - cbind(z,1-z) %*% beta, rho, sigma.sq)
  #print(cost)
  z_0 <- rep(1, each = n)
  beta_0 <- xQy(z_0, y, rho, sigma.sq)/xQy(z_0, z_0, rho, sigma.sq)
  cost_0 <- xQy(y - z_0 * beta_0, y-z_0 * beta_0, rho, sigma.sq)
  #print(cost_0)
  cost_fin <- (cost_0 - cost)
  MISE <- mean((y - cbind(z,1-z) %*% beta)^2)
  return(cost_fin)
}

rfgls_cost_t2 <- function(c_0,y,x, rho, sigma.sq){
  if(c_0 <= min(x)){
    cost_fin <- 0
  }
  if(c_0 > max(x)){
    cost_fin <- 0
  }
  if(c_0 > min(x) & c_0 <= max(x)){
    cost_fin <- rfgls_cost_t1(c_0,y,x, rho, sigma.sq)
  }
  return(cost_fin/length(y))
}

cost_fn <- function(c_0,y,x){
  if(sum(as.numeric(x > c_0)) <= 1){
    cost <- 0
  }
  if(sum(as.numeric(x <= c_0)) <= 1){
    cost <- 0
  }
  
  if(sum(as.numeric(x <= c_0)) > 1 & sum(as.numeric(x > c_0)) > 1){
    cost <- (length(y) - 1)/length(y) * var(y) - (sum(as.numeric(x > c_0))-1)/length(y) * var(y[which(x > c_0)]) -  (sum(as.numeric(x <= c_0))-1)/length(y) * var(y[which(x <= c_0)])
  }
  return(cost)
}

xQy <- function(x,y,rho, sigma.sq){
  n <- length(x)
  p1 <- x[1]*y[1] - rho * x[2] * y[1]
  pn <- x[n]*y[n] - rho * x[n-1] * y[n]
  pmid <- sum((1+rho^2) * x[2:(n-1)] * y[2:(n-1)] - rho * y[2:(n-1)] * (x[1:(n-2)] + x[3:n]))
  return((p1 + pn+ pmid)/sigma.sq)
}

rfgls_beta_t1 <- function(c_0,y,x, rho, sigma.sq){
  n <- length(y)
  z <- as.numeric(x > c_0)
  zQz <- matrix(0,2,2)
  zQz[1,1] <- xQy(z, z, rho, sigma.sq)
  zQz[1,2] <- xQy(z, 1-z, rho, sigma.sq)
  zQz[2,1] <- xQy(1-z, z, rho, sigma.sq)
  zQz[2,2] <- xQy(1-z, 1-z, rho, sigma.sq)
  zQy <- matrix(0,2,1)
  zQy[1,1] <- xQy(z, y, rho, sigma.sq)
  zQy[2,1] <- xQy(1-z, y, rho, sigma.sq)
  beta <- solve(zQz, zQy)
  return(beta)
}

rf_beta <- function(c_0,y,x){
  return(c(mean(y[which(x > c_0)]), mean(y[which(x <= c_0)])))
}

x_0 <- 0.5
n <- 1000
rho <- 0.9
sigma.sq <- 10
len <- 1:n

c_0_seq <- seq(0,1, by = 0.001)
rfgls_value_matrix <- matrix(0, length(c_0_seq), 100)
rf_value_matrix <- matrix(0, length(c_0_seq), 100)
rfgls_beta_matrix <- matrix(0,2,100)
rf_beta_matrix <- matrix(0,2,100)
rfgls_opt_vec <- rep(0, 100)
rf_opt_vec <- rep(0, 100)


Sigma <- as.matrix(rho^dist(len))+diag(n)
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}
seed_iter <- 2
set.seed(seed_iter)
eps <-sqrt(sigma.sq) * rmvn(100, rep(0,n), Sigma)
for(iter in 1:100){
  set.seed(seed_iter+iter+1)
  x <- runif(n)
  y <- as.numeric(x > x_0) + 1.5 * as.numeric(x <= x_0) + eps[,iter]
  rfgls_value_matrix[,iter] <-  sapply(c_0_seq, rfgls_cost_t2, y, x, rho, sigma.sq = 1)
  rf_value_matrix[,iter] <-  sapply(c_0_seq, cost_fn, y, x)
  rfgls_opt <- c_0_seq[which(rfgls_value_matrix[,iter] == max(rfgls_value_matrix[,iter]))]
  rfgls_beta_matrix[,iter] <- rfgls_beta_t1(rfgls_opt[1], y, x, rho, sigma.sq)
  rf_opt <- c_0_seq[which(rf_value_matrix[,iter] == max(rf_value_matrix[,iter]))]
  rfgls_opt_vec[iter] <- rfgls_opt[1]
  rf_opt_vec[iter] <- rf_opt[1]
  rf_beta_matrix[,iter] <- rf_beta(rf_opt[1], y, x)
}

  