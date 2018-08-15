d.lo=-0.5
d.hi=0.5

to.stat.1=TOSTone.raw(m=0.2, 
                      mu=0,
                      sd=1.38,
                      n=10,
                      low_eqbound=d.lo, 
                      high_eqbound=d.hi,
                      alpha=0.05,
                      plot=TRUE)

p.tost.1=max(to.stat.1$TOST_p1,to.stat.1$TOST_p2)

p.delta.2=sgpv(est.lo=to.stat.1$LL_CI_TTEST, est.hi=to.stat.1$UL_CI_TTEST, 
               null.lo=d.lo, null.hi=d.hi)
p.tost.1
p.delta.2

(to.stat.1$UL_CI_TTEST - to.stat.1$LL_CI_TTEST) / (d.hi - d.lo)