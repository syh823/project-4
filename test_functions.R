rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <-  function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- - 4*k*th[1]
  h
}


newt(c(3,5), rb, gb, hb, k = 3)

x <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240)

nll <- function(theta,x,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i)) 
  ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * x) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood 
  } ## nll

gll <- function(theta,x,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model 
  alpha <- theta[1];beta <- theta[2] ## enhances readability 
  ebt <- exp(beta*x) ## avoid computing twice 
  -c(sum(y)/alpha - sum(ebt), ## -dl/dalpha
  sum(y*x) - alpha*sum(x*ebt)) ## -dl/dbeta 
  } ## gll

hll <- function(theta,x,y) {
  ## Hessian of -ve log lik of Poisson AIDS early epidemic model 
  alpha <- theta[1];beta <- theta[2] ## enhances readability 
  ebt <- exp(beta*x) ## avoid computing twice
  H <- matrix(0,2,2) ## matrix for Hessian of -ve ll
  H[1,1] <- sum(y)/alpha^2
  H[2,2] <- alpha*sum(x^2*ebt)
  H[1,2] <- H[2,1] <- sum(x*ebt)
  H
} ## hll

newt(c(10, 10), nll, gll, hll, x = x, y = y, maxit = 200)
newt(c(2, 20), nll, gll, hll, x = x, y = y)
