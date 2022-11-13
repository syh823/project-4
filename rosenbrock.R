# Use of deriv function to create the derivatives and the Hessian matrix
grb <- deriv(expression(k * (z - x^2)^2 + (1 - x)^2),c("z", "x"), ## dif wrt these
             function(z,x,k){}, hessian = TRUE)

# Rosenbrock function where hessian attribute exists
rosenbrock1 <- function(theta, k){
  # Split input vector to z and x for deriv() usage
  z <- theta[1]; x <- theta[2]
  f <- k * (z - x^2)^2 + (1 - x)^2
  # Assigning gradient and hessian attributes to f
  attr(f, 'gradient') <- c(attr(grb(z,x,k), "gradient"))
  attr(f, 'hessian') <- matrix(attr(grb(z,x,k), "hessian"),nrow = length(theta), ncol = length(theta))
  return(f)
}

# Rosenbrock function where hessian attribute does not exist
rosenbrock2 <- function(theta, k){
  # Split input vector to z and x for deriv() usage
  z <- theta[1]; x <- theta[2]
  f <- k * (z - x^2)^2 + (1 - x)^2
  # Assigning gradient and hessian attributes to f
  attr(f, 'gradient') <- c(attr(grb(z,x,k), "gradient"))
  return(f)
}

