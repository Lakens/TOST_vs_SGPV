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

pTOSTtwo_converter=function(p_tost, 
                         m1,m2, 
                         sd1,sd2, 
                         n1,n2, 
                         alpha, 
                         low_eqbound, 
                         high_eqbound, 
                         bound,
                         var.equal){

  # the standard error for the CI won't be computed the same way 
  # if the homoscedasticity assumption is met or not 

    if(var.equal=="TRUE"){
    pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    sd <- pooled_sd*sqrt(1/n1+1/n2)
  } else {sd<-sqrt(sd1^1/n1+sd2^2/n2)}
  
  #------------------------------------------------------------------------
  # Computing the SGPV when the CI around I is non-point CI (ie lb != ub)
  #------------------------------------------------------------------------

    if (sd!=0){  # data are variable (>< constant)

      ### Computing the CI around (m1-m2) 
   
           
      # If var.equal=TRUE

       
      if (var.equal==TRUE){
        
       # computing t, based on p_tost and df
        t <- abs((qt(p_tost, df = n1+n2-2)))
        dist <- t*pooled_sd*sqrt(1/n1+1/n2) # how far is I from the closest eqbound?

       # computing I
          if (bound=="low") {
            if (p_tost <= 0.5){I = low_eqbound + dist}
            if (p_tost > 0.5){I = low_eqbound - dist}} 
          if (bound=="high") {
            if (p_tost <= 0.5){I = high_eqbound - dist}
            if (p_tost > 0.5){I = high_eqbound + dist}}
        
       # Computing the CI around I
        lb <- I - qt(1 - alpha/2, df = n1+n2-2) * pooled_sd*sqrt(1/n1+1/n2) # lower bound
        ub <- I + qt(1 - alpha/2, df = n1+n2-2) * pooled_sd*sqrt(1/n1+1/n2) # upper bound

        
      # If var.equal=FALSE
        

       } else { # if var.equal=FALSE
          t <- abs((qt(p_tost, df = df_w)))
          df_w<-(sd1^2/n1+sd2^2/n2)^2/(sd1^4/(n1^2*(n1-1))+sd2^4/(n2^2*(n2-1)))
          dist <- t*sqrt(sd1^2/n1+sd2^2/n2) # how far is I from the closest eqbound?
          
        # computing I
        if (bound=="low") {
          if (p_tost <= 0.5){I = low_eqbound + dist}
          if (p_tost > 0.5){I = low_eqbound - dist}} 
        if (bound=="high") {
          if (p_tost <= 0.5){I = high_eqbound - dist} 
          if (p_tost > 0.5){I = high_eqbound + dist}}
          
        # Computing the CI around I
        lb <- I - qt(1 - alpha/2, df = df_w) * sqrt(sd1^2/n1+sd2^2/n2) # lower bound
        ub <- I + qt(1 - alpha/2, df = df_w) * sqrt(sd1^2/n1+sd2^2/n2) # upper bound
        }

        
      ### Computing SGPV 
 
           
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
        } else {SGPV=NULL
        print("ERROR: p_tost can not be < .5 with point H0")}}    # with point H0 p1+p2=1. Because p_tost is always (max (p1,p2), no value smaller than .5 is possible)
      
      return(SGPV)
      
  }     

  #------------------------------------------------------------------------
  # Computing the SGPV when the CI around I is point CI (ie lb == ub)
  #------------------------------------------------------------------------
    if (sd==0){    # the data are constant, not variable
    # similar result, Whatever H0 is point or non-point H0
    
    print ("p_tost can not be computed because t = inf")

    # computing I
    I = m1-m2

    if (I <= high_eqbound & I >= low_eqbound){SGPV = 1  # I is included inside H0 (when non-point H0) OR I = H0 (when point H0)
    } else {SGPV = 0}                                    # I is out of H0 (when non-point H0)/ I != H0 (when point H0)
    print (paste("based on the eqb
                 ounds and I, SGPV =",SGPV))
    }
     
}
 
############ TESTS

m1<-147
m2 <- 146
sd1 <- 2
sd2<-2
n1 <- 10
n2<-10
low_eqbound = 1 
high_eqbound = 1 
alpha = 0.05
var.equal=TRUE

res<-TOSTtwo.raw(m1=m1,m2=m2,sd1=sd1,sd2=sd2,n1=n1,n2=n2,low_eqbound=low_eqbound,high_eqbound=high_eqbound,alpha=alpha,var.equal=TRUE)

SGPV1=p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)

if(res$TOST_p1>res$TOST_p2){bound="low"
}else {bound="high"}

SGPV2=pTOSTtwo_converter(p_tost=max(res$TOST_p1,res$TOST_p2),
                      bound=bound,
                      sd1=sd1,sd2=sd2,
                      n1=n1,n2=n2,
                      alpha=alpha,
                      low_eqbound=low_eqbound,
                      high_eqbound=high_eqbound,
                      var.equal=TRUE)

SGPV1
SGPV2

