options(scipen = 99)

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

# Recreate examples in Figure 2 just the check that functions work as they should ----
#1
res <- TOSTone.raw(m = 146, 
                   mu = 146,
                   sd = 500, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)

#2
res <- TOSTone.raw(m = 145.5, 
                   mu = 146,
                   sd = 250, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)

#3
res <- TOSTone.raw(m = 145, 
                   mu = 146,
                   sd = 1250, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)


#4
res <- TOSTone.raw(m = 146, 
                   mu = 146,
                   sd = 2250, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)

#5
res <- TOSTone.raw(m = 144, 
                   mu = 146,
                   sd = 1000, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)

#6
m <- 143.5
mu <- 146
sd <- 500
n <- 1000000
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha
)

p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)

#7
m <- 142
mu <- 146
sd <- 1000
n <- 1000000
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha
)

p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)

#8
m <- 141
mu <- 146
sd <- 500
n <- 1000000
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha
)

p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)



#Example where TOST + NHST is more informative than sgpv (both TOST as NHST have high p-values, sgpv is uninformatively 1)----
m <- 145
mu <- 146
sd <- 500
n <- 1000000
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha
)

p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)



step = 0.01

p_tost_list <- numeric(length(seq(140, 146, step)))
sgpv_list <- numeric(length(seq(140, 146, step)))
p_list <- numeric(length(seq(140, 146, step)))
t_list <- numeric(length(seq(140, 146, step)))

count <- 0

for(i in seq(140, 146, step)){
  count <- count + 1
  m <- i
  mu <- 146
  sd <- 0.8
  n <- 100
  low_eqbound = -3 
  high_eqbound = 3 
  alpha = 0.05
   
  invisible(capture.output(res <- TOSTone.raw(m = m, 
                                              mu = mu,
                                              sd = sd, 
                                              n = n, 
                                              low_eqbound = low_eqbound, 
                                              high_eqbound = high_eqbound, 
                                              alpha = alpha,
                                              plot = FALSE
  )))
  t <- (m - mu)/(sd/sqrt(n))
  t_list[count] <- t
  sgpv_list[count] <- p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)
  p_tost_list[count] <- max(res$TOST_p1, res$TOST_p2)
  p_list[count] <- 2 * pt(-abs(t), df = n-1)
}

plot(sgpv_list, p_tost_list)
plot(sgpv_list, p_tost_list, xlim=c(0,0.05), ylim=c(0.95, 1))
plot(sgpv_list, p_tost_list, xlim=c(0.95,1), ylim=c(0.00, 0.05))


plot(sgpv_list, p_list)
plot(p_tost_list, p_list)

sgpv_list * 2/t_list


sgpv_list-(1-p_tost_list)
plot(sgpv_list-(1-p_tost_list))

plot(sgpv_list, type="l", col = "blue")
lines(1-p_tost_list)

qnorm(p_tost_list)


#3
res <- TOSTone.raw(m = 145, 
                   mu = 146,
                   sd = 1250, 
                   n = 1000000, 
                   low_eqbound = -2, 
                   high_eqbound = 2, 
                   alpha = 0.05
)

p_delta(146+res$LL_CI_TTEST, 146+res$UL_CI_TTEST, 144, 148)

lb <- 142
ub <- 148
delta_lb <- 144
delta_ub <- 148

(min(ub, delta_ub) - max(lb, delta_lb)) / min(ub - lb, 2 * (delta_ub - delta_lb))

# equivalence bound = 4, 2/3 of the CI overlaps, SGPV is 0.666

lb <- 143
ub <- 149
delta_lb <- 144
delta_ub <- 148

(min(ub, delta_ub) - max(lb, delta_lb)) / min(ub - lb, 2 * (delta_ub - delta_lb))

# still .666! Because SGPV doesn't care you are wrong on both ends. 

lb <- 142
ub <- 145
delta_lb <- 144
delta_ub <- 148

(min(ub, delta_ub) - max(lb, delta_lb)) / min(ub - lb, 2 * (delta_ub - delta_lb))

# still .666! Because SGPV doesn't care you are wrong on both ends. 

m <- 144.5
mu <- 146
sd <- 500
n <- 1000000
low_eqbound = -2 
high_eqbound = 2 
alpha = 0.05

res <- TOSTone.raw(m = m, 
                   mu = mu,
                   sd = sd, 
                   n = n, 
                   low_eqbound = low_eqbound, 
                   high_eqbound = high_eqbound, 
                   alpha = alpha
)

p_delta(mu+res$LL_CI_TTEST, mu+res$UL_CI_TTEST, mu+low_eqbound, mu+high_eqbound)

t <- (m - mu - low_eqbound)/(sd/sqrt(n))
p <- pt(t, n-1, lower.tail = FALSE)

#TOST
(m - mu)-(1.64*(sd/sqrt(n)))
(m - mu)+(1.64*(sd/sqrt(n)))

#SGPV
ll_ci <- (m - mu)-(1.96*(sd/sqrt(n)))
ul_ci <- (m - mu)+(1.96*(sd/sqrt(n)))

(min(ul_ci, high_eqbound) - max(ll_ci, low_eqbound)) / min(ul_ci - ll_ci, 2 * (high_eqbound - low_eqbound))




(((m - mu)+(1.96*(sd/sqrt(n)))) - low_eqbound) / (((m - mu)+(1.96*(sd/sqrt(n)))) - ((m - mu)-(1.96*(sd/sqrt(n)))))

(6-3)/(6-2)
(6-3)/(6-2)


6/(6-2)-3/(6-2)

6/(6-2)

1/4 * 6

1/(6-2) * 6


1.5/ (1/6 / 1/2) 

6/12 * 6/2


  
6/4
1/4 * 6

6/5
1/5 * 6
6/14 * 7/2



(2*6) / (6*6) 
3/6 * 6/2


6/6/2


 (3)/(2)



#Calculate SGPV
lb <- mu+res$LL_CI_TTEST
ub <- mu+res$UL_CI_TTEST
delta_lb <- 144
delta_ub <- 148
(min(ub, delta_ub) - max(lb, delta_lb)) / min(ub - lb, 2 * (delta_ub - delta_lb))
p_delta(res$LL_CI_TTEST, res$UL_CI_TTEST, low_eqbound, high_eqbound)
