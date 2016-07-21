I was asked to do a Gene Set Enrichment Analysis (GSEA) for RNA-seq data.
One of the most popular tool is [GSEA](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html) from broad Institute. To better
understand the underlying method of `GSEA`, I read the [original paper: Gene set enrichment analysis: A knowledge-based
approach for interpreting genome-wide expression profiles](http://software.broadinstitute.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf) and searched on [biostars](https://www.biostars.org/p/132575/).

I quote from the biostar post:  
>so, to run GSEA you have your list of genes (L) and two conditions (or more), i.e. a microarray with normal and tumor samples. the first thing that GSEA does is to rank the genes in L based on "how well they divide the conditions" using the probe intensity values. at this point you have a list L ranked from 1...n.  
now you want to see whether the genes present in a gene set (S) are at the top or at the bottom of your list...or if they are just spread around randomly. to do that GSEA calculates the famous enrichment score, that becomes normalized enrichment score (NES) when correcting for multiple testing (FDR).  
a positive NES will indicate that genes in set S will be mostly represented at the top of your list L. a negative NES will indicate that the genes in the set S will be mostly at the bottom of your list L.  
let's say that S1 has positive NES and S2 has negative NES. let's say also that your list of 1000 genes is ordered form the most upregulated (top: 1,2,3,....) to the most downregulated (bottom: ....n-3,n-2,n-1,n). a positive NES for S1 will mean that genes over-represented in that gene set are upregulated in your dataset. negative NES for S2 instead indicated the opposite.  
in the results you will also find a heatmap the subset of you data that belong to the signature analyzed. generally what I saw is that the more significantly enriched is the gene set, the better the division between the two conditions in the heatmap.

One can run GSEA in two modes:

#### using raw gene expression data
Supply a expression data file in various [formats](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
one of the common format is GCT: Gene Cluster Text file format (*.gct)).
![](https://cloud.githubusercontent.com/assets/4106146/16968303/8f65e8e6-4dd3-11e6-9a98-093eb0bd1e86.png) 
and a phenotype label file :
![](https://cloud.githubusercontent.com/assets/4106146/16968346/ca846ed4-4dd3-11e6-89a7-be32c0c62e3b.png)

Then GSEA will calculate the rank of the genes by different [matrics](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking).

I read the [mannual](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Run_GSEA_Page) of GSEA and found:
>Metric for ranking genes. GSEA ranks the genes in the expression dataset and then analyzes that ranked list of genes. 
Use this parameter to select the metric used to score and rank the genes; use the Gene list sorting mode parameter to determine 
whether to sort the genes using the real (default) or absolute value of the metric score; and use the Gene list ordering mode 
parameter to determine whether to sort the genes in descending (default) or ascending order. 
For descriptions of the ranking metrics, see Metrics for Ranking Genes.

>Note: The default metric for ranking genes is the signal-to-noise ratio. 
To use this metric, your phenotype file must define at least two categorical phenotypes and your expression dataset must contain at least **three (3)** samples for each phenotype. If you are using a continuous phenotype or your expression dataset contains fewer than three samples per phenotype, you must choose a different ranking metric. 

>**If your expression dataset contains only one sample, you must rank the genes and use the GSEAPreranked Page to analyze the ranked list; none of the GSEA metrics for ranking genes can be used to rank genes based on a single sample.**

#### [using a pre-ranked gene list](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_GSEAPreranked_Page)



