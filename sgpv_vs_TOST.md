Equivalence Testing and the Second Generation P-Value
================
Daniël Lakens & Marie Delacre
1 juli 2018

The second generation p-value (SGPV) is the proportion of data-supported hypotheses that are also null hypotheses (Blume, McGowan, Dupont, & Greevy, (2018). The authors note that: "Using second-generation p-values can only improve rigor, reproducibility and transparency across science." It was quickly noted on Twitter that the SGPV is similar to equivalence testing (<https://twitter.com/statsepi/status/997759878503550976>, <https://twitter.com/lakens/status/995171827692515328>).

In the plot below *p*-values are calculated for the TOST equivalence testing procedure where a true population mean ranging from 140 to 150 is compared to the test value of 145 in a one-sample equivalence test where equivalence bounds are set to difference of -2 and +2 around the test value of 145. In other words, the equivalence range in the test contains all means between 143 and 147. FOr clarity, sample sizes consist of 1000000 observations. The population standard deviation is set to 500.

![](sgpv_vs_TOST_files/figure-markdown_github/sgpv_tost-1.png)

For ease of comparison we can reverse the SGPV (by calculating 1-SGPV). We see that

![](sgpv_vs_TOST_files/figure-markdown_github/1-sgpv_tost-1.png)

When the population mean is 145, and we are testing against equivalence bounds of 143 and 147 using the TOST procedure for a one-sample *t*-test with a sample size of 1000000 and a standard deviation of 500, the equivalence test is significant, *t*(999999) = 4, *p* = 0.0000317. Because the 95% CI falls completely within the equivalence bounds, the SGPV is 1.

One the other hand, if the observed mean is 140, the equivalence test is not surprisingly not significant (the observed mean is far outside the equivalence range of 143 to 147), , *t*(999999) = -6, *p* = 1. Because the 95% CI falls completely within the equivalence bounds, the SGPV is 0.

SGPV as a uniform measure of overlap
------------------------------------

It is clear the SGPV and the p-value from TOST are closely related. We can think of the SGPV as a straight line that will always overlap the p-value from an equivalence test in 3 points. When the TOST p-value is 0.5, the SGPV is also 0.5. The SGPV is 50% when the observed mean falls exactly on the lower or upper equivalence bound. When the observed mean equals the equivalence bound, the difference between the mean in the data and the equivalence bound is 0, the *t*-value for the equivalence test is also so, and thus the *p*-value is 0.5.

![](sgpv_vs_TOST_files/figure-markdown_github/unnamed-chunk-6-1.png)

Two other points always have to overlap. When the 95% CI falls completely, but only just inside the equivalence region, the TOST (which relies on a one-sided test) should be significant at an alpha level of 0.025. When the SGPV changes from 0.9999 to exactly 1 the 95% CI just touches the equivalence bound (see situation B in the plot above, where the 95% CI falls completely inside the equivalence bounds) the TOST *p*-value is 0.025. The third point where the SGPV and the *p*-value from the TOST procedure should overlap is where the SGPV changes from a positive value (i.e., 0.0001) to 0 (when the 95% CI completely falls outside of the equivalence bound, see situation C in the plot above). When the 95% CI touches the outside of the equivalence bound and the TOST *p*-value will be 0.975.

The confidence interval width is a uniformly distributed across the mean differences, in the sense that as the true mean in a one-sample t-test gets closer to the test value (in the plot below, from situation A to D, the mean gets closer to the test value by 0.1) the difference in the overlap is stable.

![](sgpv_vs_TOST_files/figure-markdown_github/unnamed-chunk-8-1.png)

For example, the SGPV from A to D is 0.7551064, 0.8061277, 0.857149, and 0.9081703. The difference in the percentage of overlap between A and B (-0.0510213) is identical to the difference in the percentage of overlap between C and D as the mean gets 0.1 closer to the test value (-0.0510213).

    ## [1] 0.07616203

As we move the means closer to the test value in steps of 0.1 across A to D the p-value calculated for normally distributed data is not uniformly distributed. The probability of observing data more extreme than the upper bound of 2 is (from A to D) 0.4207403, 0.3445783, 0.0807567, and 0.0547993. As we can see, the difference between A and B (0.076162) is not the same as the difference between C And D (0.0259574). Indeed, the difference in *p*-values is the largest as you start at *p* = 0.5 (when the observed mean falls on the test value), which is why the line in Figure 1 is the steepest at p = 0.5. Note that where the SGPV reaches 1 or 0, *p*-values closely approximate 0 and 1, but never reach these values.

References
==========

Blume, J. D., McGowan, L. D., Dupont, W. D., & Greevy, R. A. (2018). Second-generation p-values: Improved rigor, reproducibility, & transparency in statistical analyses. PLOS ONE, 13(3), e0188299. <https://doi.org/10.1371/journal.pone.0188299>
