#######

install.packages("TOSTER")
library("TOSTER")

n=40
breaks=c(6,10,20)
sims=1000
d.lo=-0.5*0.75
d.hi=0.5*0.75

keep=matrix(-99,nrow=sims,ncol=6)

for (i in 1:sims) {

z=rnorm(n,mean=0,sd=1)

to.stat.1=TOSTone.raw(m=mean(z[1:breaks[1]]), mu=0, 
	sd=sd(z[1:breaks[1]]), n=length(z[1:breaks[1]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot="FALSE")

to.stat.2=TOSTone.raw(m=mean(z[1:breaks[2]]), mu=0, 
	sd=sd(z[1:breaks[2]]), n=length(z[1:breaks[2]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot="FALSE")

to.stat.3=TOSTone.raw(m=mean(z[1:breaks[3]]), mu=0, 
	sd=sd(z[1:breaks[3]]), n=length(z[1:breaks[3]]), 
	low_eqbound=d.lo, high_eqbound=d.hi, alpha=0.05,plot="FALSE")

p.tost.1=max(to.stat.1$TOST_p1,to.stat.1$TOST_p2)
p.tost.2=max(to.stat.2$TOST_p1,to.stat.2$TOST_p2)
p.tost.3=max(to.stat.3$TOST_p1,to.stat.3$TOST_p2)

tt.stat.1=t.test(z[1:breaks[1]])$conf.int
tt.stat.2=t.test(z[1:breaks[2]])$conf.int
tt.stat.3=t.test(z[1:breaks[3]])$conf.int
	
p.delta.1=sgpv(est.lo=tt.stat.1[1], est.hi=tt.stat.1[2], 
			null.lo=d.lo, null.hi=d.hi, 
			na.constant=1e-5, length.warn=TRUE)

p.delta.2=sgpv(est.lo=tt.stat.2[1], est.hi=tt.stat.2[2], 
			null.lo=d.lo, null.hi=d.hi, 
			na.constant=1e-5, length.warn=TRUE)

p.delta.3=sgpv(est.lo=tt.stat.3[1], est.hi=tt.stat.3[2], 
			null.lo=d.lo, null.hi=d.hi, 
			na.constant=1e-5, length.warn=TRUE)

keep[i,]=cbind(p.delta.1,p.tost.1,p.delta.2,p.tost.2,p.delta.3,p.tost.3)
		}

plot(keep[,1],keep[,2],xlim=c(0,1),ylim=c(0,1),type="n",
		xlab="SGPV",ylab="TOST PV")

abline(0,1,lty=2,col="grey")

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