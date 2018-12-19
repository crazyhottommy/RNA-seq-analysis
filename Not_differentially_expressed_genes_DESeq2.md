**UPDATE: 2018-12-19**

[Mike Love](https://github.com/crazyhottommy/RNA-seq-analysis/issues/8) pointed out that in DESeq2,

>Around 2014 we implemented `lfcThreshold` as an argument of `results()` which takes care of everything for the use

Please refer to [this section of the vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tests-of-log2-fold-change-above-or-below-a-threshold) for the "not DE" analysis.


I read a [post](http://watson.nci.nih.gov/~sdavis/blog/testing-for-non-differentially-expressed-genes/) by Sean Davis @NCI discussing
how to detect non-differential expressed genes.
I keep it a note here as it is sometimes desirable to find such genes or ChIP-seq regions that are not changed after treatment.


Simon Anders answered the question quite thoroughly in text like so. And I quote:
>One of the new features of DESeq2 is that it provides a empirical-Bayes style shrinkage estimates of coefficients (i.e., log fold changes) and also estimates standard errors for these coefficients (taken from the the reciprocal curvature of the posterior). Your task is one of the application we had in mind for this. You want to know whether the fold change is zero or nearly zero, and you also want to be ascertained that this estimate of a nearly-zero fold change is a precise one. So, to find genes whose absolute log fold change is _reliably_ smaller than some threshold, take all those genes for which the estimated abs log fold change is below the threshold and the standard error is well below the threshold, too. The reason, why I write "below a threshold" rather than "equal to zero" is that it's biologically unreasonable to assume that a treatment has exactly zero effect on a gene's expression. Gene regulation is such an interconnected network that, at least in my view, every gene will react ever so slightly to every perturbance, but often, this reaction is too small to be measurable. This is why you should define a threshold for "biologically unlikely to be relevant". Let's say, we don't believe that an expression change of less than 15% (0.2 on a log2 scale) is worth bothering. So, if you might take all genes with abs log2 fold change below 0.2, and furthermore require that the standard error of this log2 fold change estimate is below, say, 0.05, than you will get gene with true values well below at least 0.3 or so. In the end, the thresholds won't matter too much but if you want real p values, this is possible, too: Simply divide the log2 fold change estimates by their standard errors to get z scores (this works because the sampling distribution of our fold change estimates is reasonably close to normal), and do two one-sided tests (TOST).

> To recap, you want to test the null hypothesis that your log fold change beta is larger than some threshold theta, i.e., you want to find genes for which you can reject this null hypothesis and hence say with some certainty that |beta| < theta. The TOST scheme suggests, as its name implies to perform two pone sided tests, one with the null hypothesis H0_A: beta > theta, the other with H0_B: beta < -theta, and then use the larger of the two p values.

```{r}
# Let's use the example data from the pasilla package
library( DESeq2 )
library( pasilla )
data( "pasillaGenes" )

# Create a DESeq2 data object from the pasilla data
dse <- DESeqSummarizedExperimentFromMatrix(
   counts(pasillaGenes), pData(pasillaGenes)[,c("condition","type")],
   ~ type + condition )

# Perform a standard DESeq2 analysis
dse <- DESeq(dse)

# The log2 fold changes are found here
beta <- results(dse)$log2FoldChange

# Just to make the following clearer, I should point out that the
# "results" accessor is just a short-cut for this access to the rowData:
all( beta == mcols(rowData(dse))$conditionuntreated, na.rm=TRUE )
# (returns TRUE)

# The log fold change estimates all come with standard error information
# which we find in the rowData (maybe we should copy this to the
# 'results', too)
betaSE <- mcols(rowData(dse))$SE_conditionuntreated

# Internally, the Wald test is implemented as a simple two-sided
# z test of beta/betaSE. Two demonstrate this, we to the test
# manually and compare
pvalDE <- 2 * pnorm( abs( beta ), sd = betaSE, lower.tail=FALSE )
all( abs( pvalDE - results(dse)$pvalue ) < 1e-15, na.rm=TRUE )
# (returns TRUE)

# This was the test for DE, of course, i.e., small pvalDE means that
# the gene's expression change (the true value of beta) is not zero

# What we want is the opposite, namely find gene, for which abs(beta)
# is smaller than some threshold, theta
theta <- .3

# So, we do our two one-sided tests. For a one-sided z test, we
# simply use tail probabilities from the normal distribution.

# First, the test of H0_A: true_beta > thr
pA <- pnorm( beta, thr, betaSE, lower.tail=TRUE )

# Next, the test of H0_B: true_beta < -thr
pB <- pnorm( beta, -thr, betaSE, lower.tail=FALSE )

# The overall p value is the maximum, because we want to reject H0_A
# and H0_B simultaneously
pvalTOST <- pmax( pA, pB )


# Let's adjust our two p values with BH:
sigDE <- p.adjust( pvalDE, "BH" ) < .1
sigSmall <- p.adjust( pvalTOST, "BH" ) < .1

# And make an MA plot, with sigDE in red and sigSmall in green
plot(
   rowMeans( counts(dse,normalized=TRUE) ), beta,
   log="x", pch=20, cex=.2,
   col = 1 + sigDE + 2*sigSmall )
# Plot is attached.

```


