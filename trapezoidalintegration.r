##### Trapezoidal Integration
##### By Eric Cai - The Chemical Statistician
##### define the function for trapezoidal integration

trapezoidal.integration = function(x, f)
{
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }
  
  # check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand "f" is not numeric.')
  }
  
  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(f))
  {
    stop('The lengths of the variable of integration and the integrand do not match.')
  }
  
  ### finish checks
  
  # obtain length of variable of integration and integrand
  n = length(x)
  
  # integrate using the trapezoidal rule
  integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))
  
  # print the definite integral
  return(integral)
}