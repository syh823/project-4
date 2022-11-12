# separate function for finite difference

newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6){

  # original = func(theta)
  # grad (theta)
  # if is.nan(original) == T or is.nan(grad(theta)) == T: stop(the objective or derivatives are not finite at the initial theta)
  # maxit = 0
  # while grad not absolute value less than tol times 
      # the absolute value of the objective function plus fscale or maxit <= 100
    # maxit = maxit + 1
    # if hess not provided, find hess
      # use finite difference page 72 notes
    # check if hess if positive definite using chol
      # if not 
        # keep adding identity matrix until chol can be computed (use try)
    # delta = -(hess(theta))^(-1)*grad(theta)
    # inverse of hess use chol and backsolve
    # newtheta = theta + delta
    # new = func(newtheta)
    # max.half = 0
      # while new > original or max.half <= 20
      # max.half = max.half + 1
      # newtheta = theta + delta/2
      # new = func(newtheta)
    # if max.half > 20 stop("step fails to reduce the objective despite trying max.half step halvings")
    # grad(newtheta)
  # if maxit > 100 stop()
  # if hess(theta) not posdef warning("hess matrix not posdef")
}







