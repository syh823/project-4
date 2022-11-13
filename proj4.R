# separate function for finite difference

newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6){

  # original = func(theta)
  original=func(theta)
  
  # grad (theta)
  origin_grad=grad(theta)
  
  # if is.nan(original) == T or is.nan(grad(theta)) == T: stop(the objective or derivatives are not finite at the initial theta)
  if(is.nan(original) == T | sum(is.nan(origin_grad))!=0)
    stop("the objective or derivatives are not finite at the initial theta")

  # if hess not provided, find hess
  if(is.null(hess)==TRUE){
    # use finite difference page 72 notes
    hess=XXXX
  }

  # maxit = 0
  maxit=0
  
  # while grad not absolute value less than tol times 
  while(sum(abs(grad(theta))>tol)!=0 & maxit <= 100){ # the absolute value of the objective function plus fscale 
    # maxit = maxit + 1
    maxit=maxit+1
    
    # check if hess if positive definite using chol
    flag=try (chol(hess(theta)))

    # if not 
    if(class(flag)[1]=='try-error'){
      # keep adding identity matrix until chol can be computed (use try)
    }
    
    # inverse of hess use chol and backsolve
    inverse_hess=solve(hess(theta))
    
    # delta = -(hess(theta))^(-1)*grad(theta)
    delta = -inverse_hess%*%grad(theta)
    
    # newtheta = theta + delta
    newtheta=theta+delta
    
    # new = func(newtheta)
    new = func(newtheta)
    
    # max.half = 0
    max.half=0
    
    # while new > original and max.half <= 20
    while(new>original& max.half<=20){
      # max.half = max.half + 1
      max.half=max.half+1
      
      # newtheta = theta + delta/2
      newtheta= theta+ delta/2
      # new = func(newtheta)
      
      new=func(newtheta)
    }
    # if max.half > 20 stop("step fails to reduce the objective despite trying max.half step halvings")
    if(max.half > 20)
      stop("step fails to reduce the objective despite trying max.half step halvings")
    # grad(newtheta)
    theta=newtheta
  }
  # if maxit > 100 stop()
  if(maxit>100)
    stop("fail")

  # if hess(theta) not posdef warning("hess matrix not posdef")
  if(class(try (chol(hess(theta))))[1]=='try-error')
    warning("hess matrix not posdef")
  print(theta)
  print(maxit)
}

newt(c(-2.3,-3.9),rb,gb,hb)





