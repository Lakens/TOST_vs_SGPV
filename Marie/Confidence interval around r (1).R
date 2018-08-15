for (package in c("psych","parallel","smoothmest","fGarch","TOSTER")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

#--------------------------------------------------------------------------

# The distribution of r is skewed when r!=0
# when r = 0, the sampling distribution or r is symmetric
# when r gets closer of +-1, the sampling distribution of r is increasingly skewed.
# For this reason, a correction for the bias is required.
# Commonly, reasearchers use the Fisher's z transformation in order to compute CI:
  # 1) convert r into z: z=log((1+r)/(1-r))/2.    z has a std error of 1/sqrt(n-3)
  # 2) computing z bounds: bound = z+-qnorm(1-alpha/2)*1/sqrt(n-3)  
  # 3) convert z bounds into r bounds: rbound= (exp(2*z)-1)/(1+exp(2*z))

#---------------------------------------------------------------------------------
 
# The following section shows that r distribution is skewed, while z distribution is normal

correlDV <- function(n, correlation=.5,nb_DV=4,sd) {
  # generate nb_DV variables, stored in a matrix
  X=matrix(0,n,nb_DV)  
  for(v in 1:nb_DV){
    X[,v] <- rnorm(n=n,m=0,sd=sd)
  }
  # generate the expected correlation matrix    
  C <- matrix(correlation, nrow = nb_DV, ncol = nb_DV)
  diag(C) <- 1 
  L <- chol(C) # Choleski factorization
  # Return a lower triangular matrix Z such as C=t(Z)%*%Z   
  # induce correlation, without changing X1 (as long as Z is a lower triangular matrix)
  DV <- X %*% L
  return(DV)
}

samplingcor_plot=function(nSims,n,correlation,sd){
  
  count=0
  cor_list<-rep(0,length(1:nSims))
  z_list<-rep(0,length(1:nSims))
  
  proportion_r<-matrix(0,length(seq(0,.9,.1)),1)
  colnames(proportion_r)=paste("proportion, when n=",n)
  rownames(proportion_r)=paste("r >=",seq(0,.9,.1))
  
  for (i in 1:nSims){
    count=count+1
    
    data=correlDV(n=n, 
                  correlation=correlation,
                  nb_DV=2,
                  sd=sd) 
    cor_list[count]=cor(data[,1],data[,2])  
    z_list[count]=log((1+cor_list[count])/(1-cor_list[count]))/2
  }
  
  par(xpd=T,mar=c(4,2,2,2))
  plot(density(cor_list),main="",xlab=paste("std error or r=",round(sd(cor_list),3),"\n","std error of z=",round(sd(z_list),3)),xlim=c(-3,3))
  lines(density(z_list),lty=2,col="blue")
  legend("topright",pch=c(15,15),lty=c(1,2),col=c("black","blue"),legend=c("r distribution","z distribution"))
  abline(v=1,lty=2,col="grey")
  abline(v=-1,lty=2,col="grey")
  
  row=0
  for (j in seq(0,.9,.1)){
    row=row+1
    proportion_r[row,]=round(sum(abs(cor_list)>=j)/length(cor_list),3)}
  print(proportion_r)
  
}

# Illustration

samplingcor_plot(nSims=100000,n=30,correlation=0,sd=2) # not skewed
samplingcor_plot(nSims=100000,n=100,correlation=.5,sd=2) # not skewed
samplingcor_plot(nSims=100000,n=100,correlation=.2,sd=2) # not skewed
samplingcor_plot(nSims=100000,n=100,correlation=.8,sd=2) # not skewed

#---------------------------------------------------------------------------------

# The following section show that after the Fisher's z transformation, the bounds
# of the CI around correlation will always fall inside the interval [-1;1].
# There is therefore no need to truncate the CI, as with proportions tests

n=10
alpha=.05
step=.01

binf_list=rep(0,length(seq(-1,1,step)))
bsup_list=rep(0,length(seq(-1,1,step)))
binf_z_list=rep(0,length(seq(-1,1,step)))
bsup_z_list=rep(0,length(seq(-1,1,step)))
i_list=rep(0,length(seq(-1,1,step)))
count=0

for (i in seq(-.99,.99,step)){
  count=count+1
  # convert r into z
  r=i
  z=log((1+r)/(1-r))/2
  
  # computing z bounds
  binf_z_list[count]=z-qnorm(1-alpha/2)/sqrt(n-3)
  bsup_z_list[count]=z+qnorm(1-alpha/2)/sqrt(n-3)

  # convert z bounds into r bounds
  binf_list[count]=(exp(2*binf_z_list[count])-1)/(1+exp(2*binf_z_list[count]))
  bsup_list[count]=(exp(2*bsup_z_list[count])-1)/(1+exp(2*bsup_z_list[count]))

  i_list[count]=i  
}

min(binf_z_list)
max(binf_z_list)
 
min(bsup_z_list)
max(bsup_z_list)
