# Input parameters

# In order to be able to convert p_tost into SGPV, there are many parameters we need to know:

# From which eqbound is computed p_tost?
	# if p_tost = TOST_p1 (= p_tost is computed from the low_eqbound), then bound="low" 
      # else if p_tost = TOST_p2 (= p-tost is computed from the high_eqbound), then bound="high" 
# n, sd 
# alpha
# low_eqbound and high_eqbound

#-----------------------------------------------------------------------------------------------

pTOST_converter=function(p_tost, 
                         sd, 
                         n, 
                         alpha, 
                         low_eqbound, 
                         high_eqbound, 
                         bound){
  
  # computing t, based on p_tost
  t <- abs((qt(p_tost, df = n-1))) 
  
  # computing I
  dist <- t * sd/sqrt(n) # how far is I from the closest eqbound?
  
  #LAKENS: I I fixed this now - you can't know based on non-sig if bound is left or right from alpha. 
  #It should be based on whether p is larger or smaller than 0.5. I smaller than 0.5, bound - dist. 
  if (bound=="low") {
    if (p_tost < 0.5){I = low_eqbound + dist}
    if (p_tost > 0.5){I = low_eqbound - dist}} 
  
  if (bound=="high") {
    if (p_tost < 0.5){I = high_eqbound - dist}
    if (p_tost > 0.5){I = high_eqbound + dist}}
  
  # Computing the CI around I
  lb <- I - qt(1 - alpha/2, df = n - 1) * sd/sqrt(n) # lower bound
  ub <- I + qt(1 - alpha/2, df = n - 1) * sd/sqrt(n) # upper bound
  
  if (lb > high_eqbound | ub < low_eqbound){SGPV = 0} 
  else if(lb > low_eqbound & ub < high_eqbound){SGPV = 1} 
  else {SGPV <- (min(ub, high_eqbound) - max(lb, low_eqbound))/min(ub - lb, 2 * (high_eqbound - low_eqbound))}
  
  return(SGPV)
}




m <- 143.9
mu <- 146
m-mu
sd <- 3
n <- 36
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m,mu = mu,sd = sd,n = n,low_eqbound = low_eqbound,high_eqbound = high_eqbound,alpha = alpha)
pTOST_converter(p_tost=res$TOST_p1,bound="low",sd=3,n=36,alpha=.05,low_eqbound=-2,high_eqbound = 2)

pTOST_converter(p_tost=res$TOST_p2,
                bound = "high",
                sd = 3,
                n = 36,
                alpha = .05,
                low_eqbound = -2,
                high_eqbound = 2)

