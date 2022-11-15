#   Newton's method of optimization to find the minimum of a given function 
# starts with applying a second order Taylor Series approximation at an initial
# guess of theta, the vector of parameters to be optimized. This approximation 
# reguires the gradient vector and square Hessian matrix, the first and second  
# order partial derivatives of the given function with respect to theta 
# respectively. The next guess of theta is then the minimum of the quadratic 
# approximation at the 
# previous guess, checking that the Hessian matrix at that guess is positive 
# definite. If it isn't, it is perturbed by adding a multiple of the identity 
# matrix, large enough to force positive definiteness. To check that the new  
# guess is getting closer to the minimum of the function, the function value at 
# the new guess must be smaller than the value at the previous guess. If this 
# condition is not met, the difference between the new guess and the previous 
# one is halved until it is met. 

# This algorithm is repeated until convergence is reached, which is when the 
# gradient at theta is approximately 0 and its Hessian matrix is 
# positive definite.


# separate function for finite difference
Hfd = function(theta, grad, eps, ...){
  grad0 <- grad(theta,...)
  Hfd <- matrix(0,length(theta),length(theta)) ## approximated Hessian matrix
  for (i in 1:length(theta)) { 
    ## finite difference approximation for each parameter (by row)
    th1 <- theta; th1[i] <- th1[i] + eps ## increase theta by small increment eps 
    grad1 <- grad(th1,...) ## compute resulting grad at th1
    Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivatives by row
  }
  return(Hfd)
}


newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=200,eps=1e-6){
  
  # original = func(theta)
  original_f = func(theta, ...)
  
  # grad (theta)
  grad_vec = grad(theta, ...)
  
  # if initial guess is the minimum
  if(is.null(hess)==TRUE) hess_mat = Hfd(theta = theta, grad = grad, eps = eps, ...)
  else hess_mat = hess(theta,...)
  
  # if is.nan(original) == T or is.nan(grad(theta)) == T: stop(the objective or derivatives are not finite at the initial theta)
  if(is.nan(original_f) == T | sum(is.nan(grad_vec))!=0)
    stop("the objective or derivatives are not finite at the initial theta")
  
  
  
  # maxit = 0
  count_maxit=0
  
  # while grad not absolute value less than tol times 
  while(sum(abs(grad_vec)>tol*abs(original_f+fscale))!=0 & count_maxit < maxit){ # the absolute value of the objective function plus fscale 
    # maxit = maxit + 1
    count_maxit=count_maxit+1
    # if hess not provided, find hess
    if(is.null(hess)==TRUE) hess_mat = Hfd(theta = theta, grad = grad, eps = eps, ...)
    else hess_mat = hess(theta,...)
    # check if hess if positive definite using chol
    flag=try (chol(hess_mat), silent = T)
    
    # if not 
    if(class(flag)[1]=='try-error'){
      # keep adding identity matrix until chol can be computed (use try)
      #d = -min(eigen(hess_mat)$values)+0.1
      #hess_mat = hess_mat + diag(d, nrow(hess_mat))
      # hess_mat=nearPD(hess_mat)$mat
      d = 1e-6*norm(hess_mat)
      hess_mat_perturbed = hess_mat + diag(d, nrow(hess_mat))
      flag_perturbed = try (chol(hess_mat_perturbed), silent = T)
      while (class(flag_perturbed)[1]=='try-error') {
        d = d * 10
        hess_mat_perturbed = hess_mat + diag(d, nrow(hess_mat))
        flag_perturbed =try (chol(hess_mat_perturbed), silent = T)
      }
      hess_mat = hess_mat_perturbed
    }
    
    # inverse of hess use chol and backsolve
    #inverse_hess = solve(hess_mat)
    
    # delta = -(hess(theta))^(-1)*grad(theta)
    #delta=-inverse_hess%*%grad(theta,...)
    delta = backsolve(chol(hess_mat),forwardsolve(t(chol(hess_mat)),-grad_vec))
    
    # newtheta = theta + delta
    newtheta=theta+delta
    
    # new = func(newtheta)
    new_f = func(newtheta,...)
    
    # max.half = 0
    count_max.half=0
    
    # while new > original and max.half <= 20
    while((new_f > original_f ||is.na(new_f))  && count_max.half<=max.half){
      delta = delta/2
      # max.half = max.half + 1
      count_max.half=count_max.half+1
      # print(func(newtheta,...))
      
      # newtheta = theta + delta/2
      newtheta= theta+ delta
      # new = func(newtheta)
      
      new_f = func(newtheta,...)
    }
    # if max.half > 20 stop("step fails to reduce the objective despite trying max.half step halvings")
    if(count_max.half > max.half)
      stop("step fails to reduce the objective despite trying max.half step halvings")
    # grad(newtheta)
    theta = newtheta
    original_f = new_f
    grad_vec = grad(theta,...)
  }
  
  # if maxit > 100 stop()
  if(count_maxit == maxit)
    warning("failed to achieve convergence at maxit iterations")
  
  output = list(func(theta,...),theta, count_maxit, grad_vec)
  names(output) = c("f", "theta", "iter", "g")
  
  # if hess(theta) not posdef warning("hess matrix not posdef")
  if(class(try (chol(hess_mat)))[1]=='try-error'){
    return(output)
    warning("hess matrix not posdef")
  } else {
    inverse_hess = backsolve(chol(hess_mat),forwardsolve(t(chol(hess_mat)),diag(length(theta))))
    output[[5]] = inverse_hess
    names(output) = c("f", "theta", "iter", "g", "Hi")
    return(output)
  }
  
  
 
  #output = list(func(theta,...),theta, count_maxit, grad_vec, inverse_hess)
  #names(output) = c("f", "theta", "iter", "g", "Hi")
  
  return(output)
}










