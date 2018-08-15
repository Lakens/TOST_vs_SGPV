library(TOSTER)

#Create SGPV funtion from https://github.com/LucyMcGowan/sgpvalue/blob/master/R/p_delta.R
#' Second Generation P-value

p_delta <- function(lb, ub, delta_lb, delta_ub) {
  
  # special case: infinite CI and H0 bounds in the same direction
  if ((delta_lb == -Inf & lb == -Inf) | (delta_ub == Inf & ub == Inf)) {
    return(1)
  }
  
  # usual case: non-point CI & non-point Ho
  # pdelta = |CI intersect Ho| / min{ |CI|, 2|Ho| }
  if (delta_lb != delta_ub & lb != ub) {
    if (lb > delta_ub | ub < delta_lb) {
      return(0)
    } else if(lb > delta_lb & ub < delta_ub){
      return(1)
    } else {
      return(
        (min(ub, delta_ub) - max(lb, delta_lb)) /
          min(ub - lb, 2 * (delta_ub - delta_lb))
      )
    }
  }
  # special case 1: point CI, w/ or w/out a point H0
  # pdelta = 0 if CI is inside the Ho
  # pdelta = 1 if CI is inside the Ho
  if (lb == ub) {
    if (lb <= delta_ub & lb >= delta_lb){
      return(1)
    } else {
      return(0)
    }
  }
  # special case 2: point H0 & non-point CI
  # pdelta = 1/2 if H0 is inside the CI
  # pdelta = 0 if H0 is outside the CI
  if (delta_lb == delta_ub & lb != ub) {
    if (delta_lb <= ub & delta_lb >= lb) {
      return(1/2)
    } else {
      return(0)
    }
  }
}

###Figure 1, 5 and 6: SGPV vs . TOST

# Figure 1
#n=120
#r_pop=.0
#low_eqbound_r=r_pop-.2
#high_eqbound_r=r_pop+.2
#alpha=.05
#step = 0.001

# Figure 5
#n=120
#r_pop=.0
#low_eqbound_r=r_pop-.1
#high_eqbound_r=r_pop+.1
#alpha=.05
step = 0.001

# Figure 6
#n=120
#r_pop=.5
#low_eqbound_r=r_pop-.2
#high_eqbound_r=r_pop+.2
#alpha=.05
#step = 0.001

# Figure 7
#n=120
#r_pop=.5
#low_eqbound_r=r_pop-.05
#high_eqbound_r=r_pop+.05
#alpha=.05
#step = 0.001

# Figure 8 
n=120
r_pop=.5
low_eqbound_r=r_pop-.1
high_eqbound_r=r_pop+.1
alpha=.05
step = 0.001

#Try out Daniel (I want the deviation at 0.5 to be greater)
n=30
r_pop=.0
low_eqbound_r=r_pop-.45
high_eqbound_r=r_pop+.45
alpha=.05
step = 0.01


p_tost_list <- numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
sgpv_list <- numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
p_list <- numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
t_list <- numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
i_list <- numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
count <- 0

lowbound_list<-numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))
highbound_list<-numeric(length(seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)))

for(i in seq(r_pop-2*(high_eqbound_r-r_pop),r_pop+2*(high_eqbound_r-r_pop), step)){
  count <- count + 1
  r <- i

  invisible(capture.output(res<-TOSTr(n=n,
                                      r=r,
                                      low_eqbound_r=low_eqbound_r,
                                      high_eqbound_r=high_eqbound_r,
                                      alpha=alpha,
                                      plot=FALSE)))
 
  t_list[count]=r*sqrt(n-2)/sqrt(1-r^2)
  p_list[count]=2*(1-pt(abs(t_list[count]),df=n-2))
  sgpv_list[count] <- p_delta(res$LL_CI_TTEST,res$UL_CI_TTEST, low_eqbound_r, high_eqbound_r)
  p_tost_list[count] <- max(res$TOST_p1, res$TOST_p2)
  i_list[count]=i
  lowbound_list[count]=res$LL_CI_TTEST
  highbound_list[count]=res$UL_CI_TTEST
} 

par(xpd=T)  
plot(i_list,1-sgpv_list, type="l", col = "blue",ylim=c(0,1),xlab="Correlation",ylab="",xlim=c(-1,1),cex=1.5)
lines(i_list,p_tost_list,lty=2)
abline(h=.975,lty=3,col="lightgrey")
abline(h=.025,lty=3,col="lightgrey")
abline(h=.5,lty=3,col="lightgrey")
legend(.3,1.25,lty=c(1,3),legend=c("1-sgpv","p TOST"),col=c("blue","black"),box.lty=0)
sum(sgpv_list==1)/length(sgpv_list)
sum(p_tost_list < .05)/length(p_tost_list)


#sgpv_list[28] # value=0
#p_tost_list[28] # value > .975
#sgpv_list[29] # value!=0
#p_tost_list[29] # value < .975

#sgpv_list[683] # value=0
#p_tost_list[683] # value > .975
#sgpv_list[682] # value!=0
#p_tost_list[682] # value < .975

#sgpv_list[496] # value=1
#p_tost_list[496] # value < .025
#sgpv_list[497] # value!=1
#p_tost_list[497] # value > .025

#sgpv_list[404] # value=1
#p_tost_list[404] # value < .025
#sgpv_list[405] # value!=1
#p_tost_list[405] # value > .025



###Figure 2: 

n=30
low_eqbound_r=-.45
high_eqbound_r=.45
r1=-.45
r2=0
r3=.45

res1<-TOSTr(n=n,r=r1,low_eqbound_r=low_eqbound_r,high_eqbound_r=high_eqbound_r,alpha=alpha,plot=FALSE, verbose = FALSE)
res2<-TOSTr(n=n,r=r2,low_eqbound_r=low_eqbound_r,high_eqbound_r=high_eqbound_r,alpha=alpha,plot=FALSE, verbose = FALSE)
res3<-TOSTr(n=n,r=r3,low_eqbound_r=low_eqbound_r,high_eqbound_r=high_eqbound_r,alpha=alpha,plot=FALSE, verbose = FALSE)

LL95_r1=res1$LL_CI_TTEST
UL95_r1=res1$UL_CI_TTEST
LL95_r2=res2$LL_CI_TTEST
UL95_r2=res2$UL_CI_TTEST
LL95_r3=res3$LL_CI_TTEST
UL95_r3=res3$UL_CI_TTEST

par(xpd=F)
plot(NA, ylim=c(0,.7), xlim=c(-1,1), bty="l", yaxt="n", ylab="",xlab="Correlation")
abline(v = c(r1, r2, r3),
       lty = 2,
       col = "grey")

points(x=r1, y=0.6, pch=15, cex=1.1)
segments(LL95_r1,0.6,UL95_r1,0.6, lwd=2)
points(x=r2, y=0.4, pch=15, cex=1.1)
segments(LL95_r2,0.4,UL95_r2,0.4, lwd=2)
points(x=r3, y=0.2, pch=15, cex=1.1)
segments(LL95_r3,0.2,UL95_r3,0.2, lwd=2)

###Figure 3: 

par(xpd=F)
r_list=rep(0,length(seq(-.99,.99,.01)))
z_list=rep(0,length(seq(-.99,.99,.01)))
count=0

for (i in seq(-.99,.99,.01)){
  count=count+1
  r_list[count]=i
  z_list[count]=log((1+i)/(1-i))/2
  }

plot(r_list,z_list,xlab="correlation r", ylab="Fisher's z score",pch=20,cex=.5)
x=seq(-.99,.99,.01)
lines(x,x,col="red")
abline(h=0)

###Figure 4:

r_list=numeric(length(seq(-3.5,3.5,.01)))
z_list=numeric(length(seq(-3.5,3.5,.01)))
count=0


for (i in seq(-3.5,3.5,.01)){
  count=count+1
  z_list[count]=i
  r_list[count]=(exp(2*z_list[count])-1)/(1+exp(2*z_list[count]))
}

plot(z_list,r_list,xlab="Fisher's z score", ylab="correlation r",pch=1,cex.pch=.5)
abline(h=0)
x=seq(-3.5,3.5,.01)
lines(x,x,col="red")
 