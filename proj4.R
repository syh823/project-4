# separate function for finite difference
# separate function for finite difference
Hfd = function(theta, grad, eps, ...){
  grad0 <- grad(theta,...)
  Hfd <- matrix(0,length(theta),length(theta)) ## finite diference Hessian 
  for (i in 1:length(theta)) { ## loop over parameters
    th1 <- theta; th1[i] <- th1[i] + eps ## increase theta by eps 
    grad1 <- grad(th1,...) ## compute resulting grad at th1
    Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivs by row
  }
  return(Hfd)
}


newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6){

  # original = func(theta)
  original = func(theta, ...)
  
  # grad (theta)
  grad_vec = grad(theta, ...)
  
  # if is.nan(original) == T or is.nan(grad(theta)) == T: stop(the objective or derivatives are not finite at the initial theta)
  if(is.nan(original) == T | sum(is.nan(grad_vec))!=0)
    stop("the objective or derivatives are not finite at the initial theta")



  # maxit = 0
  maxit=0
  
  # while grad not absolute value less than tol times 
  while(sum(abs(grad_vec)>tol)!=0 & maxit <= 100){ # the absolute value of the objective function plus fscale 
    # maxit = maxit + 1
    maxit=maxit+1
    # if hess not provided, find hess
    if(is.null(hess)==TRUE) hess_mat = Hfd(theta = theta, grad = grad, eps = eps, ...)
    else hess_mat = hess(theta)
    # check if hess if positive definite using chol
    flag=try (chol(hess_mat))

    # if not 
    if(class(flag)[1]=='try-error'){
      # keep adding identity matrix until chol can be computed (use try)
    }
    
    # inverse of hess use chol and backsolve
    inverse_hess = solve(hess_mat)
    
    # delta = -(hess(theta))^(-1)*grad(theta)
    delta = -inverse_hess%*%grad_vec
    
    # newtheta = theta + delta
    newtheta=theta+delta
    
    # new = func(newtheta)
    new = func(newtheta)
    
    # max.half = 0
    max.half=0
    
    # while new > original and max.half <= 20
    while(new > original& max.half<=20){
      # max.half = max.half + 1
      max.half=max.half+1
      
      # newtheta = theta + delta/2
      newtheta= theta+ delta/2
      # new = func(newtheta)
      
      new = func(newtheta)
    }
    # if max.half > 20 stop("step fails to reduce the objective despite trying max.half step halvings")
    if(max.half > 20)
      stop("step fails to reduce the objective despite trying max.half step halvings")
    # grad(newtheta)
    theta = newtheta
    grad_vec = grad(theta)
  }
  # if maxit > 100 stop()
  if(maxit>100)
    stop("fail")

  # if hess(theta) not posdef warning("hess matrix not posdef")
  if(class(try (chol(hess_mat)))[1]=='try-error')
    warning("hess matrix not posdef")
  print(theta)
  print(maxit)
}







