# Input parameters

# In order to be able to convert p_tost into SGPV, there are many parameters we need to know:

# From which eqbound is computed p_tost?
# if p_tost = TOST_p1 (= p_tost is computed from the low_eqbound), then bound="low" 
# else if p_tost = TOST_p2 (= p-tost is computed from the high_eqbound), then bound="high" 
# n, sd 
# alpha
# low_eqbound and high_eqbound

#-----------------------------------------------------------------------------------------------
library(TOSTER)

pTOSTone_converter=function(p_tost, 
                         sd, 
                         n, 
                         alpha, 
                         low_eqbound, 
                         high_eqbound, 
                         bound){

  #------------------------------------------------------------------------
  # Computing the SGPV when the CI around I is non-point CI (ie lb != ub)
  #------------------------------------------------------------------------
    if (sd!=0){  # data are variable (>< constant)

      # computing t, based on p_tost
      t <- abs((qt(p_tost, df = n-1))) 
      
      # computing I
      dist <- t * sd/sqrt(n) # how far is I from the closest eqbound?
      if (bound=="low") {
        if (p_tost <= 0.5){I = low_eqbound + dist}
        if (p_tost > 0.5){I = low_eqbound - dist}} 
      
      if (bound=="high") {
        if (p_tost <= 0.5){I = high_eqbound - dist}
        if (p_tost > 0.5){I = high_eqbound + dist}}
      
      # Computing the CI around I
      lb <- I - qt(1 - alpha/2, df = n - 1) * sd/sqrt(n) # lower bound
      ub <- I + qt(1 - alpha/2, df = n - 1) * sd/sqrt(n) # upper bound
      
      qt(1 - alpha/2, df = 8.20001736626575)

      # Computing the SGPV
      
      # If H0 is non-point H0
      if(low_eqbound!=high_eqbound){ 
        if (p_tost>(1-alpha/2)){SGPV = 0} 
        else if(p_tost<(alpha/2)){SGPV = 1} 
        else {SGPV <- (min(ub, high_eqbound) - max(lb, low_eqbound))/min(ub - lb, 2 * (high_eqbound - low_eqbound))}} 
      
      # If H0 is point H0
      else {  
        if (p_tost >=.5 & p_tost <= 1-alpha/2){SGPV = 1/2}          # i.e. lb <= low_eqbound & low_eqbound <= ub
                                                                    # point H0 is inside the CI around I
        else if (p_tost>1-alpha/2){SGPV = 0                         # point H0 is outside the CI around I
        }else {SGPV=NULL
          print("ERROR: p_tost can not be < .5 with point H0")}}    # with point H0 p1+p2=1. Because p_tost is always (max (p1,p2), no value smaller than .5 is possible)
      
      return(SGPV)
      }     

  #------------------------------------------------------------------------
  # Computing the SGPV when the CI around I is point CI (ie lb == ub)
  #------------------------------------------------------------------------
    if (sd==0 ){    # the data are constant, not variable
    # similar result, Whatever H0 is point or non-point H0
    
    print ("p_tost can not be computed because t = inf")

    # computing I
    I = m-mu

    if (I <= high_eqbound & I >= low_eqbound){SGPV = 1  # I is included inside H0 (when non-point H0) OR I = H0 (when point H0)
    } else {SGPV = 0}}                                    # I is out of H0 (when non-point H0)/ I != H0 (when point H0)
    print (paste("based on the eqbounds and I, SGPV =",SGPV))    

}
 
############ TESTS

m<-146
mu <- 146
sd <- 3
n <- 36
low_eqbound = 0 
high_eqbound = 0 
alpha = 0.01

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha)


pTOSTone_converter(p_tost=.3,bound=bound,sd=sd,n=n,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound)

SGPV1=p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)

if(res$TOST_p1>res$TOST_p2){bound="low"
}else {bound="high"}

SGPV2=pTOSTone_converter(p_tost=max(res$TOST_p1,res$TOST_p2),bound=bound,sd=sd,n=n,alpha=alpha,low_eqbound=low_eqbound,high_eqbound=high_eqbound)

SGPV1
SGPV2
