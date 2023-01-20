<h1 align="center">Newton's Optimisation Method</h1> 

## Project Overview

Newton's method of optimization to find the minimum of a given function starts with applying a second order Taylor Series approximation at an initial guess of theta, the vector of parameters to be optimized. This approximation requires the gradient vector and square Hessian matrix, the first and second order partial derivatives of the given function with respect to theta respectively. The next guess of theta is then the minimum of the quadratic approximation at the  previous guess, checking that the Hessian matrix at that guess is positive definite. If it isn't, it is perturbed by adding a multiple of the identity matrix, large enough to force positive definiteness. To check that the new  guess is getting closer to the minimum of the function, the function value at the new guess must be smaller than the value at the previous guess. If this  condition is not met, the difference between the new guess and the previous one is halved until it is met. 

This algorithm is repeated until convergence is reached, which is when the gradient at theta is approximately 0 and its Hessian matrix is positive definite. 

## Project Content

This project consists of two functions, `Hfd` and `newt`. `Hfd` is a separate function applies finite difference approximations to derive Hessian matrix to be used in `newt`. `newt` is designed such that the output is comparable to the built-in `nlm` function. The documentation of arguments for each function is as follows:

- `Hfd(theta, grad, eps, ...)`
  - theta: parameters of the objective function
  - grad: gradient function of the objective function
  - eps: a properly small number to differentiate gradient (for approximation)
  - ...: passed arguments from objective function

- `newt(theta, func, grad, hess = NULL, ..., tol = 1e-8, fscale = 1, maxit = 100, max.half = 200, eps = 1e-6)`
  - theta: an initial guess of optimization parameters
  - func: the objective function to minimize
  - grad: the gradient function of objective function
  - hess: the Hessian matrix of objective function. It may be provided or not.
  - ...: passed arguments from objective function
  - tol,fscale: values used to test convergence
  - maxit: maximum number of iterations
  - max.half: maximum number of times of halving a step before giving up at that step 
  - eps: a properly small number to be used in finite difference formulas

The final marked 15/18 (83%) project can be found [here](https://github.com/syh823/project-4/blob/main/proj4.R)
The revised version after comments can be found [here]()
