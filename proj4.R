# separate function for finite difference

newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6){
  # start loop
  # find func value at theta = original
  # while grad not absolute value less than tol times 
    # the absolute value of the objective function plus fscale or maxit == 100
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
    # while new > original or max.half == 20
    # max.half = max.half + 1
    # newtheta = theta + delta/2
  # 
}