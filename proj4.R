# 1. Nurfahimah binti Mohd Ghazali, s2464388@ed.ac.uk
# 2. Yanming Gu, s2304572@ed.ac.uk
# 3. Yuheng Song, s2447118@ed.ac.uk

# Github address: https://github.com/syh823/project-4

# All three of us worked together on the project. Contributions are equal for
# all group members.


#   Newton's method of optimization to find the minimum of a given function 
# starts with applying a second order Taylor Series approximation at an initial
# guess of theta, the vector of parameters to be optimized. This approximation 
# requires the gradient vector and square Hessian matrix, the first and second  
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



# Separate function: using finite difference approximations to derive Hessian 
# matrix by differencing the gradient passed from newt function
# Called when Hessian matrix is not supplied in newt function.
Hfd = function(theta, grad, eps, ...){
  # theta: parameters of the function
  # grad: gradient of the function
  # eps: a properly small number to differentiate gradient 
  grad0 <- grad(theta,...)
  if (length(theta)==1){
    Hfd = matrix((grad(theta + eps,...) - grad0)/eps, 1, 1)
  } else{
    Hfd <- matrix(0,length(theta),length(theta)) ## approximated Hessian matrix
    for (i in 1:length(theta)) { 
      ## finite difference approximation for each parameter (by row)
      th1 <- theta; th1[i] <- th1[i] + eps ## increase theta by small increment eps 
      grad1 <- grad(th1,...) ## compute resulting grad at th1
      Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivatives by row
    }
  }
  return(Hfd)
}



# This is the main function of Newton Optimizer
# Given a function and its gradient with an initial guess of theta, we can find
# the minimization of the objective function within limited trials.
newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=200,eps=1e-6){
  # The explaination of every each input:
    # theta: an initial guess of optimization parameters
    # func: the objective function to minimize
    # grad: the gradient of objective function
    # hess: the Hessian matrix of objective function. It may be provided or not.
    # tol,fscale: both them are used to test convergence
    # maxit: maximum number of iterations
    # max.half: maximum number of times of halving a step before giving up that step 
    # eps: a properly small number to be used in finite difference formulas
  
  original_f = func(theta, ...)  # Objective function at initial value of theta
  grad_vec = grad(theta, ...)  # Gradient at initial value of theta
  
  # If Hessian matrix is not provided, then we need to call our function to calculate it
  if(is.null(hess)==TRUE) hess_mat = Hfd(theta = theta, grad = grad, eps = eps, ...)
  else hess_mat = hess(theta,...)
  
  # Check the objective or derivatives at initial theta are both not finite
  if(is.nan(original_f) == T | sum(is.nan(grad_vec))!=0)
    stop("the objective or derivatives are not finite at the initial theta")
  
  # Count the number of iterations we've tried
  count_maxit=0
  
  # Check the condition on convergence
  while(sum(abs(grad_vec)>tol*abs(original_f+fscale))!=0 & count_maxit < maxit){ 
    count_maxit=count_maxit+1 ## One more iteration
    # If Hessian matrix is not provided, calculate it
    if(is.null(hess)==TRUE) hess_mat = Hfd(theta = theta, grad = grad, eps = eps, ...)
    else hess_mat = hess(theta,...)
    # Check if Hessian matrix is positive definite using Cholesky decomposition
    flag=try (chol(hess_mat), silent = T)
    
    # If Hessian matrix is not positive definite, we perturb it by add a multiple 
    # of the identity matrix
    if(class(flag)[1]=='try-error'){
      d = 1e-6*norm(hess_mat)
      hess_mat_perturbed = hess_mat + diag(d, nrow(hess_mat))
      flag_perturbed = try (chol(hess_mat_perturbed), silent = T)
      # Keep perturb Hessian matrix until it can be Cholesky decomposed
      while (class(flag_perturbed)[1]=='try-error') {
        d = d * 10
        hess_mat_perturbed = hess_mat + diag(d, nrow(hess_mat))
        flag_perturbed =try (chol(hess_mat_perturbed), silent = T)
      }
      hess_mat = hess_mat_perturbed
    }
    
    # Calculate the step towards minimizing the objective function
    delta = backsolve(chol(hess_mat),forwardsolve(t(chol(hess_mat)),-grad_vec))
    
    # Update theta and objective function
    newtheta=theta+delta
    new_f = func(newtheta,...)
    
    # Count the number of a step being halved
    count_max.half=0
    
    # Check if the step optimize the objective function within max.half times
    while((new_f > original_f ||is.finite(new_f) == FALSE)  && count_max.half<=max.half){
      delta = delta/2  # Halve the step is it can not reduce the function
      count_max.half=count_max.half+1 
      # Update theta and objective function
      newtheta= theta+ delta  
      new_f = func(newtheta,...)
    }
    # The step fails to reduce the objective despite trying max.half step halvings
    if(count_max.half > max.half)
      stop("step fails to reduce the objective despite trying", as.character(max.half), "step halvings")
    
    # Update theta, objective function and its gradient
    theta = newtheta
    original_f = new_f
    grad_vec = grad(theta,...)
  }
  
  # Fail to minimize the function within maxit times of iterations
  if(count_maxit == maxit)
    warning("failed to achieve convergence at ", as.character(maxit), " iterations")
  
  # Output values at the minimum
  output = list(func(theta,...),theta, count_maxit, grad_vec)
  names(output) = c("f", "theta", "iter", "g")
  
  # Check if Hessian matrix is positive definite at convergence
  # If not, then give a warning
  if(class(try (chol(hess_mat)))[1]=='try-error'){
    return(output) 
    warning("hess matrix not posdef")  
  } else {
    # Hessian is positive definite, then the inverse Hessian could also be returned.
    inverse_hess = chol2inv(chol(hess_mat))
    output[[5]] = inverse_hess
    names(output) = c("f", "theta", "iter", "g", "Hi")
    return(output)
  }
  
}







