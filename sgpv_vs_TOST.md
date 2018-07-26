Equivalence Testing and the Second Generation P-Value
================
Daniël Lakens & Marie Delacre
26 juli 2018

The RMarkdown file with the reproducible code of this text is [here](https://github.com/Lakens/TOST_vs_SGPV/blob/master/sgpv_vs_TOST.Rmd).

The second generation *p*-value (SGPV) is a new descriptive statistic that was recently proposed to "improve rigor, reproducibility and transparency across science" (Blume, McGowan, Dupont, & Greevy, (2018). The SGPV is 'the proportion of data-supported hypotheses that are also null hypotheses'. The researcher specify an equivalence range around the null hypothesis that specifies values that are considered practically equivalent to the null-hypothesis. The SGPV is the proportion of the confidence interval (CI) around the observed effect estimate that falls within this equivalence range. If the CI falls completely inside the equivalence range the SGPV is 1, and if the CI falls completely outside the equivalence rnge the SGPV is 0.

The SGPV has strong similarities with an already existing approach known as equivalence testing (Lakens, 2017; Rogers, Howard, & Vessey, 1993). In an equivalence test, an equivalence range is based upon a smallest effect size of interest. In the Two One-Sided Tests (TOST) approach to equivalence testing data is tested against the upper and lower bounds of the equivalence range (e.g., a difference of -2 and +2). If both one-sided tests reject the presence of effects more extreme than the equivalence bounds when can act as if there is no meaningful effect. Because two one-sided tests are performed, equivalence can be declared when a 1-2\*alpha CI (e.g., when the alpha level is 0.05, a 90% CI) falls completely within the equivalence range of -2 and +2.

Surprisingly, Blume et al (2018) do not discuss equivalence testing in their article, despite the strong conceptual similarities. Here, we aim to examine the similarities and differences between equivalence testing using the TOST procedure and the SGPV. Our goal is to allow researchers to choose the statistic that best answers the question they are interested in.

The relationship between *p*values from TOST and SGPV
=====================================================

In the plot below *p*-values are calculated for the TOST equivalence testing procedure where a true population mean ranging from 140 to 150 is compared to the test value of 145 in a one-sample equivalence test where equivalence bounds are set to difference of -2 and +2 around the test value of 145. In other words, the equivalence range in the test contains all means between 143 and 147. Blume et al (2018) rely on the z-distribution, while to TOST package uses the *t*-distribution (which is more accurate at smaller sample sizes). To make sure the SGPV give basically identical results, sample sizes consist of 1000000 observations (for which the *t*-distribution and z-distribution are basically identical). The population standard deviation is set to 500 to still give some variation in responses. Our conclusions should hold to the same extend for more realistic numbers (e.g., N = 100, SD = 1).

![](sgpv_vs_TOST_files/figure-markdown_github/sgpv_tost-1.png) *Figure 1*: Comparison of *p*-values from TOST (black line) and SGPV (dotted grey line) across a range of true population means (x-axis) tested against a mean of 145 in a one-sample *t*-test with a sample size of 1000000 and a standard deviation of 500.

The SGPV treats the equivalence range as the null-hypothesis, while the TOST procedure treats the values outside of the equivalence range as the null-hypothesis. For ease of comparison we can reverse the SGPV (by calculating 1-SGPV) to make the two tests more comparable. We see that the *p*-value from the TOST procedure and the SGPV follow each other closely. ![](sgpv_vs_TOST_files/figure-markdown_github/1-sgpv_tost-1.png)

*Figure 2*: Comparison of *p*-values from TOST (black line) and 1-SGPV (dotted grey line) across a range of true population means (x-axis) tested against a mean of 145 in a one-sample *t*-test with a sample size of 1000000 and a standard deviation of 500.

When the population mean is 145 and we are testing against equivalence bounds of 143 and 147 using the TOST procedure for a one-sample *t*-test with a sample size of 1000000 and a standard deviation of 500, the equivalence test is significant, *t*(999999) = 4, *p* = 0.0000317. Because the 95% CI falls completely within the equivalence bounds, the SGPV is 1 (see Figure 1).

One the other hand, if the observed mean is 140, the equivalence test is not significant (the observed mean is far outside the equivalence range of 143 to 147), *t*(999999) = -6, *p* = 1 (or more accuratelty, *p* &gt; .999 as *p*-values are bounded between 0 and 1). Because the 95% CI falls completely outside the equivalence bounds, the SGPV is 0 (see Figure 1).

SGPV as a uniform measure of overlap
------------------------------------

It is clear the SGPV and the *p*-value from TOST are closely related. We can think of the SGPV as a straight line that will always overlap the *p*-value from an equivalence test in 3 points. When the TOST *p*-value is 0.5, the SGPV is also 0.5 (note that the reverse is not true). The SGPV is 50% when the observed mean falls exactly on the lower or upper equivalence bound. When the observed mean equals the equivalence bound, the difference between the mean in the data and the equivalence bound is 0, the *t*-value for the equivalence test is also 0, and thus the *p*-value is 0.5 (situation A).

![](sgpv_vs_TOST_files/figure-markdown_github/unnamed-chunk-6-1.png) *Figure 3*: Means, normal distribution, and 95% CI for three example datasets that illustrate the relationship between *p*-values from TOST and SGPV.

Two other points always have to overlap. When the 95% CI falls completely, but only just inside the equivalence region, the TOST (which relies on a one-sided test) should be significant at an alpha level of 0.025. When the SGPV changes from &lt;1 to 1 the 95% CI is exactly equal to one of the equivalence bounds (see situation B in the plot above, where the 95% CI falls completely inside the equivalence bounds) the TOST *p*-value is 0.025. The third point where the SGPV and the *p*-value from the TOST procedure should overlap is where the SGPV changes from a positive value (i.e., 0.0001) to 0 (when the 95% CI completely falls outside of the equivalence bound, see situation C in the plot above). When the 95% CI touches the outside of the equivalence bound and the TOST *p*-value will be 0.975.

The confidence interval width is a uniformly distributed across the mean differences, in the sense that as the true mean in a one-sample t-test gets closer to the test value (in the plot below, from situation A to D, the mean gets closer to the test value by 0.1) the difference in the overlap is stable.

![](sgpv_vs_TOST_files/figure-markdown_github/unnamed-chunk-8-1.png) *Figure 4*: Means, normal distribution, and 95% CI for data with a sample size of 1000000 and a standard deviation of 500 for samples where the true population mean is 1.5, 1.4, 1.3, and 1.2.

For example, the SGPV from A to D is 0.7551064, 0.8061277, 0.857149, and 0.9081703. The difference in the percentage of overlap between A and B (-0.0510213) is identical to the difference in the percentage of overlap between C and D as the mean gets 0.1 closer to the test value (-0.0510213).

As we move the means closer to the test value in steps of 0.1 across A to D the *p*-value calculated for normally distributed data is not uniformly distributed. The probability of observing data more extreme than the upper bound of 2 is (from A to D) 0.1586553, 0.1150697, 0.0807567, and 0.0547993. As we can see, the difference between A and B (0.0435856) is not the same as the difference between C And D (0.0259574). Indeed, the difference in *p*-values is the largest as you start at *p* = 0.5 (when the observed mean falls on the test value), which is why the line in Figure 1 is the steepest at *p* = 0.5. Note that where the SGPV reaches 1 or 0, *p*-values closely approximate 0 and 1, but never reach these values.

What are the Relative Strengths and Weaknesses of Equivalence Testing and SGPV?
-------------------------------------------------------------------------------

Given the strong relationship between SGPV and equivalence testing, a logical question is to ask what the introduction of SGPV adds to the existing statistical approaches, including equivalence tests, and what the relative strengths and weaknesses of either approach are. First of all, SGPV is a descriptive statistic (unlike the *p*-value that is calculated for an equivalence test, which is an inferential statistic). It numerically summarizes the information that is visually present in a plot (such as Figure 3) displaying the equivalence range and the 95% CI around the observed effect.

The SGPV is 1 for tests where *p*-values for the TOST procedure differ. For example, different equivalence tests with *p* = 0.024 and *p* = 0.0001 have a SGPV of 1. Although a SGPV of 1 or 0 has a clear interpretation (we can reject effects outside or inside the equivalence range) intermediate values are not as easy to interpret (e.g., it is unclear how we would interpret a SGPV of 0.56 versus 0.65). Since the SGPV is always directly related to a *p*-value from the TOST procedure, different SGPV can be interpreted in the same manner as different *p*-values. From a Fisherian viewpoint, the lower the *p*-value, the worse the fit of the data with a specific model, and analogously, the lower the SGPV the worse the fit of the data with the equivalence range. From a Neyman-Pearson approach to statistics, only the dichotomous rejection of values outside of the equivalence range (TOST *p* &lt; *α* or SGPV = 1) allows you to act as if the null-hypothesis is true while controlling our error rate at a known maximum.

\[To be extended\]

Discussion
==========

It seems Blume et al (2018) where not aware of the existence of equivalence tests, and we believe that our explanation of the similarities between the TOST procedure and the SGPV provides some useful context to interpret the contribution of second generation *p*values to the statistical toolbox. The novelty lies in its use as a descriptive statistic. The added benefit of calculating the proportion of overlap of a 95% CI with the equivalence range, and using this percentage to describe the data, remains somewhat unclear for practical puroposes. Nevertheless, our only goal is to clarify the relationship between a newly proposed statistic and the already existing TOST approach used to test for equivalence, and let researchers make an informed decision about which statistical approach provides the best answer to their question.

References
==========

Blume, J. D., McGowan, L. D., Dupont, W. D., & Greevy, R. A. (2018). Second-generation *p*-values: Improved rigor, reproducibility, & transparency in statistical analyses. PLOS ONE, 13(3), e0188299. <https://doi.org/10.1371/journal.pone.0188299> Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests, Correlations, and Meta-Analyses. Social Psychological and Personality Science, 8(4), 355–362. <https://doi.org/10.1177/1948550617697177> Rogers, J. L., Howard, K. I., & Vessey, J. T. (1993). Using significance tests to evaluate equivalence between two experimental groups. Psychological Bulletin, 113(3), 553–565. <http://dx.doi.org/10.1037/0033-2909.113.3.553>
