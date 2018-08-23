#######
## Adjustments to Jeffrey's code by Daniel
#install.packages("TOSTER")

library("TOSTER")

##########################################################################
# Create function for Second Generation P-value ----
sgpv <- function(est.lo, est.hi, null.lo, null.hi) {
  
  # special case: infinite CI and H0 bounds in the same direction
  if ((null.lo == -Inf & est.lo == -Inf) | (null.hi == Inf & est.hi == Inf)) {
    return(1)
  }
  
  # usual case: non-point CI & non-point Ho
  # pdelta = |CI intersect Ho| / min{ |CI|, 2|Ho| }
  if (null.lo != null.hi & est.lo != est.hi) {
    if (est.lo > null.hi | est.hi < null.lo) {
      return(0)
    } else if(est.lo > null.lo & est.hi < null.hi){
      return(1)
    } else {
      return(
        (min(est.hi, null.hi) - max(est.lo, null.lo)) /
          min(est.hi - est.lo, 2 * (null.hi - null.lo))
      )
    }
  }
  
  # special case 1: point CI, w/ or w/out a point H0
  # pdelta = 0 if CI is inside the Ho
  # pdelta = 1 if CI is inside the Ho
  if (est.lo == est.hi) {
    if (est.lo <= null.hi & est.lo >= null.lo){
      return(1)
    } else {
      return(0)
    }
  }
  
  # special case 2: point H0 & non-point CI
  # pdelta = 1/2 if H0 is inside the CI
  # pdelta = 0 if H0 is outside the CI
  if (null.lo == null.hi & est.lo != est.hi) {
    if (null.lo <= est.hi & null.lo >= est.lo) {
      return(1/2)
    } else {
      return(0)
    }
  }
}
################################################################################

#Example 1

n=1000
sd <- 10
breaks=c(100,500,1000)
sims=2000
d.lo=-0.5*0.75
d.hi=0.5*0.75

keep=matrix(-99,nrow=sims,ncol=6)

for (i in 1:sims) {

z=rnorm(n,mean=0,sd=sd)

to.stat.1=TOSTone.raw(m=mean(z[1:breaks[1]]), mu=0, 
	sd=sd(z[1:breaks[1]]), n=length(z[1:breaks[1]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot=F)

to.stat.2=TOSTone.raw(m=mean(z[1:breaks[2]]), mu=0, 
	sd=sd(z[1:breaks[2]]), n=length(z[1:breaks[2]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot=F)

to.stat.3=TOSTone.raw(m=mean(z[1:breaks[3]]), mu=0, 
	sd=sd(z[1:breaks[3]]), n=length(z[1:breaks[3]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot=F)

p.tost.1=max(to.stat.1$TOST_p1,to.stat.1$TOST_p2)
p.tost.2=max(to.stat.2$TOST_p1,to.stat.2$TOST_p2)
p.tost.3=max(to.stat.3$TOST_p1,to.stat.3$TOST_p2)

tt.stat.1=t.test(z[1:breaks[1]])$conf.int
tt.stat.2=t.test(z[1:breaks[2]])$conf.int
tt.stat.3=t.test(z[1:breaks[3]])$conf.int
	
p.delta.1=sgpv(est.lo=tt.stat.1[1], est.hi=tt.stat.1[2], 
			null.lo=d.lo, null.hi=d.hi)

p.delta.2=sgpv(est.lo=tt.stat.2[1], est.hi=tt.stat.2[2], 
			null.lo=d.lo, null.hi=d.hi)

p.delta.3=sgpv(est.lo=tt.stat.3[1], est.hi=tt.stat.3[2], 
			null.lo=d.lo, null.hi=d.hi)

keep[i,]=cbind(p.delta.1,p.tost.1,p.delta.2,p.tost.2,p.delta.3,p.tost.3)
		}

plot(keep[,1],keep[,2],xlim=c(0,1),ylim=c(0,1),type="n",
		xlab="SGPV",ylab="TOST PV")

abline(1,-1,lty=2,col="grey")

points(keep[,1],keep[,2],pch=20,col="dodgerblue")
points(keep[,3],keep[,4],pch=20,col="forestgreen")
points(keep[,5],keep[,6],pch=20,col="firebrick")

legend("topright",bty="n",pch=20,
		col=c("dodgerblue","forestgreen","firebrick"),
		c(paste("n = ",breaks[1],sep=""),
		  paste("n = ",breaks[2],sep=""),
		  paste("n = ",breaks[3],sep="")))

keep[keep[,1]==0,]


test.dat=rnorm(40,mean=0,sd=1)
to.stat.3=TOSTone.raw(m=mean(test.dat), mu=0, 
		sd=sd(test.dat), n=length(test.dat), 
		low_eqbound=-0.5, high_eqbound=0.5, alpha=0.05,plot="FALSE")


###
##
#




#Single example to see what is going on.

m <- 0.8
n <- 10
sd <- 1.5
d.lo <- -0.5*0.75
d.hi <- 0.5*0.75

z <- rnorm(n, mean = 0, sd = sd)

res <- TOSTone.raw(m = m, 
                   mu = 0,
                   sd = sd,
                   n = n,
                   low_eqbound = d.lo,
                   high_eqbound = d.hi,
                   alpha = 0.05,
                   plot = T)

tt.stat.1 <- c(res$LL_CI_TTEST,res$UL_CI_TTEST)

sgpv(est.lo = tt.stat.1[1], 
     est.hi = tt.stat.1[2],
     null.lo = d.lo,
     null.hi = d.hi)
