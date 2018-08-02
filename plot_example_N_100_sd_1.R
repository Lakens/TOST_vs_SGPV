library(TOSTER)
source("p_delta_function.R")

step = 0.01

p_tost_list <- numeric(length(seq(140, 150, step)))
sgpv_list <- numeric(length(seq(140, 150, step)))
p_list <- numeric(length(seq(140, 150, step)))
t_list <- numeric(length(seq(140, 150, step)))

count <- 0

for(i in seq(140, 150, step)){
  count <- count + 1
  m <- 140
  mu <- 145
  sd <- 1
  n <- 100
  low_eqbound = -2 
  high_eqbound = 2 
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

plot(NA, 
     ylim = c(0, 1), 
     xlim = c(0, 1001),
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "Mean")
axis(1, at = seq(0,1000,100), labels = seq(140,150,1), las = 1)
axis(2, at = seq(0,1,0.1), labels = seq(0,1,0.1), las = 1)

points(1-sgpv_list, type="l", col = "darkgrey", lwd = 3, lty = 3)
points(p_tost_list, lwd = 3)

