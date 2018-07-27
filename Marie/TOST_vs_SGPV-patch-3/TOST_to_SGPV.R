# Input parameters
# res = the stored results from the TOST analysis
# n is the sample size
# sd is the standard deviation
#-----------------------------------------------------------------------------------------------

TOST_to_SGPV <- function(tost_res){
  
  # computing t, based on p_tost
  t <- min(tost_res$TOST_t1, tost_res$TOST_t2) 
  
  n <- tost_res$TOST_df + 1
  sd <- (tost_res$diff-tost_res$low_eqbound)/tost_res$TOST_t1 * sqrt(n)
  
  # computing I
  dist <- t * sd/sqrt(n) # how far is I from the closest eqbound?
  
  #LAKENS: I I fixed this now - you can't know based on non-sig if bound is left or right from alpha. 
  #It should be based on whether p is larger or smaller than 0.5. I smaller than 0.5, bound - dist. 
  if (tost_res$TOST_p1 >= tost_res$TOST_p2) {
    if (tost_res$TOST_p1 <= 0.5){I = tost_res$low_eqbound + dist}
    if (tost_res$TOST_p1 >= 0.5){I = tost_res$low_eqbound - dist}} 
  
  if (tost_res$TOST_p2 >= tost_res$TOST_p1) {
    if (tost_res$TOST_p2 <= 0.5){I = tost_res$low_eqbound - dist}
    if (tost_res$TOST_p2 >= 0.5){I = tost_res$low_eqbound + dist}} 
  
  # Computing the CI around I
  lb <- I - qt(1 - tost_res$alpha/2, df = n - 1) * sd/sqrt(n) # lower bound
  ub <- I + qt(1 - tost_res$alpha/2, df = n - 1) * sd/sqrt(n) # upper bound
  
  if (lb > tost_res$high_eqbound | ub < tost_res$low_eqbound){
    SGPV = 0
    } 
  else if(lb > tost_res$low_eqbound & ub < tost_res$high_eqbound){
    SGPV = 1
    } 
  else {SGPV <- (min(ub, tost_res$high_eqbound) - max(lb, tost_res$low_eqbound))/min(ub - lb, 2 * (tost_res$high_eqbound - tost_res$low_eqbound))}
  
  return(SGPV)
}
