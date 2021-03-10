### Function to perform a Completely randomized experiment given the treatment group sizes

#' Completely Randomized Design (CRD)
#'
#' @description 
#' Generates an assignment under completely randomized design (CRD). 
#' @param data_frame A data frame corresponding to the full sample of units. 
#' @param n_treat Number of treatment groups.
#' @param treat_sizes A vector of treatment group sizes. If \code{control = TRUE}, 
#' the first element of \code{treat_sizes} should be the control group size.
#' @param control If \code{TRUE}, treatments are labeled as 0,1,...,g-1 (0 representing the control group). 
#' If \code{FALSE}, they are labeled as 1,2,...,g.
#' @return The original data frame augmented with the column of the treatment indicator.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), 
#' ``Randomized and Balanced Allocation of Units into Treatment Groups Using the Finite Selection Model for \code{R}''.
#' @examples
#' # Consider N = 12, n1 = n2 = n3 = 4.
#' df_sample = data.frame(index = 1:12, x = c(20,30,40,40,50,60,20,30,40,40,50,60))
#' # Draw a random assignment from CRD.
#' fc = crd(data_frame = df_sample, n_treat = 3, treat_sizes = c(4,4,4))
#' # Get vector of treatment assignments.
#' Z_crd = fc$Treat

crd = function(data_frame, n_treat, treat_sizes, control = FALSE)
{
  N = nrow(data_frame)
  Z = rep(0,N)
  unit.alloc = 1:N
  for(i in 1:n_treat)
  {
    ind = sample(unit.alloc, treat_sizes[i], replace = F)
    Z[ind] = i
    unit.alloc = setdiff(unit.alloc,ind)
  }
  
  if(control == TRUE)
  {
    Z = Z-1
  }
  data_frame.allocated = cbind(data_frame,Z)
  colnames(data_frame.allocated)[ncol(data_frame.allocated)] = 'Treat'
  return(data_frame.allocated)
}

