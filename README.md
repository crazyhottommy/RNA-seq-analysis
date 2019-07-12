# RNA-seq analysis

### General sequencing data analysis materials
* [Next-Gen Sequence Analysis Workshop (2015)](http://angus.readthedocs.org/en/2015/) held by [Titus Brown](http://genomecenter.ucdavis.edu/people/faculty/name/c-titus-brown/)  (now in UC Davis)
* [Fall 2015, BMMB 852: Applied Bioinformatics](http://www.personal.psu.edu/iua1/2015_fall_852/main_2015_fall_852.html) by        [Istvan Albert](http://www.personal.psu.edu/iua1/) from Penn state University. He developed the all-time popular 
   [biostars](https://www.biostars.org/)  
* Steven Turner in UVA is maitaining a list of training opportunities for [genomic data analysis](http://stephenturner.us/edu.html)
*  Jeff Leek group's recommended [genomic papers](https://github.com/jtleek/genomicspapers/)
* [awesome tutorial for NGS file format](http://binf.snipcademy.com/lessons/sequence-file-formats)
* [UVA Bioconnector Workshops](http://bioconnector.org/workshops/)
* [Explaining your errors QC fail](https://sequencing.qcfail.com/)
* [EMBL-EBI has a very comprehensive list of courses for online training](http://www.ebi.ac.uk/training/online/)

### RNA-seq specific 

*  [Introduction to RNA-seq analysis youtube video](https://www.youtube.com/watch?v=OEbjHPk20C0&feature=youtu.be&a)  
*  [RNAseq differential expression analysis – NGS2015](https://monsterbashseq.wordpress.com/2015/08/26/rnaseq-differential-expression-analysis-ngs2015/)  
*  [Kallisto and sleuth tutorial](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html) blazing fast RNA-seq analysis by Lior Patcher's lab.    [A sleuth for RNA-Seq](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)  
*  pathway analysis using [GAGE](https://github.com/ajwije/150826_pathway_analysis/blob/master/Tutorial_150827.Rmd) 
*  [Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview](http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html)
*  [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/14/6/671.full)
*  [RNA-seq tutorial wiki](https://github.com/crazyhottommy/rnaseq_tutorial) Informatics for RNA-seq: A web resource for analysis on the cloud.  
*  [RNA-seqlopedia](http://rnaseq.uoregon.edu/)  Great introduction of RNA-seq from sample preparation to data analysis
*  [RNAseq data analysis from data carpentry](https://github.com/datacarpentry/rnaseq-data-analysis)
*  [paper: Isoform prefiltering improves performance of count-based methods for analysis of differential transcript usage](http://www.genomebiology.com/2016/17/1/12)
*  [paper: A survey of best practices for RNA-seq data analysis](http://www.genomebiology.com/2016/17/1/13)
*  [paper: Reproducibility of high-throughput mRNA and small RNA sequencing across laboratories](http://www.nature.com/nbt/journal/v31/n11/full/nbt.2702.html)
*  [paper: Cross-platform normalization of microarray and RNA-seq data for machine learning applications](https://peerj.com/articles/1621/#results|discussion). [Tool](https://github.com/greenelab/TDM) 
*  [review: Translating RNA sequencing into clinical diagnostics: opportunities and challenges](http://www.nature.com/nrg/journal/v17/n5/full/nrg.2016.10.html)
* [paper: Thousands of large-scale RNA sequencing experiments yield a comprehensive new human gene list and reveal extensive transcriptional noise](https://www.biorxiv.org/content/early/2018/05/29/332825)

### RNA-seq experimental design 
* [Thinking about Designing RNA Seq Experiments to Measure Differential Gene Expression: The Basics](http://michelebusby.tumblr.com/post/26913184737/thinking-about-designing-rna-seq-experiments-to) a blog post 
* [Tutorial: Rna Seq Experimental Design For Measuring Differential Gene Expression from biostars](https://www.biostars.org/p/65824/)
* [Scotty - Power Analysis for RNA Seq Experiments](http://bioinformatics.bc.edu/marthlab/scotty/scotty.php)  
* [Experimental Design in Differential Abundance analysis web server](http://edda.gis.a-star.edu.sg/)  
* [Experimental Design in Differential Abundance analysis bioconductor package](http://www.bioconductor.org/packages/devel/bioc/html/EDDA.html)


### Quality Control

* [QoRTs](http://hartleys.github.io/QoRTs/): a comprehensive toolset for quality control and data processing of RNA-Seq experiments  
* [QUaCRS](http://bioserv.mps.ohio-state.edu/QuaCRS/index.php/pages/view/downloads)    
* [RSeQC](http://rseqc.sourceforge.net/) RNA-seq data QC
* [RNA-SeqQC](https://www.broadinstitute.org/cancer/cga/rna-seqc)

### Normalization, quantification, and differential expression

Normalization is essential for RNAseq analysis. However, one needs to understand the underlining assumptions for each methods. Most methods assume there is no global changes between conditions (e.g. TMM normalization). However, this may not be true when global effect occurs. For example, if you delete a gene that controls transcription, you expect to see global gene expression reduction. In that case, other normalization methods need to be considered. (e.g. spike-in controls). The same principle applies to other high-throughput sequencing data such as ChIPseq.

read this very important paper by Rafael A Irizarry: [Genome-wide repressive capacity of promoter DNA methylation is revealed through epigenomic manipulation](https://www.biorxiv.org/content/early/2018/08/01/381145)

DESseq2 normalization by Simon Anders:
>To estimate the library size, simply taking the total number of (mapped or unmapped) reads is, in our experience, not a good idea.
Sometimes, a few very strongly expressed genes are differentially expressed, and as they make up a good part of the total counts, they skew this number. After you divide by total counts, these few strongly expressed genes become equal, and the whole rest looks differentially expressed.

>The following simple alternative works much better:

>- Construct a "reference sample" by taking, for each gene, the geometric mean of the counts in all samples.

>- To get the sequencing depth of a sample relative to the reference, calculate for each gene the quotient of the counts in your sample divided by the counts of the reference sample. Now you have, for each gene, an estimate of the depth ratio. 

>- Simply take the median of all the quotients to get the relative depth of the library.

>This is what the `estimateSizeFactors` function of our DESeq package doese.

If one wants to use a set of genes that are not affected by the global change, do

```r
dds = newCountDataSet(CountTable, Design$condition )
dds <- estimateSizeFactors(dds, 
                           controlGenes = rownames(dds) %in% norm_genes)
dds_global <- estimateSizeFactors(dds)
dds_global <- DESeq(dds_global)
res_global <- results(dds_global)
```
or give self-defined size factors.

```r
sizeFactors(dds) = c(my_Values)
```
*  [A Comparison of Methods: Normalizing High-Throughput RNA Sequencing Data](http://biorxiv.org/content/early/2015/09/03/026062)
*  [Errors in RNA-Seq quantification affect genes of relevance to human disease](http://www.genomebiology.com/2015/16/1/177)  
*  [A comprehensive evaluation of ensembl, RefSeq, and UCSC annotations in the context of RNA-seq read mapping and gene quantification](http://www.biomedcentral.com/1471-2164/16/97)  
*  [Comparing the normalization methods for the differential analysis of Illumina high-throughput RNA-Seq data](http://www.biomedcentral.com/1471-2105/16/347)
*  paper: [Union Exon Based Approach for RNA-Seq Gene Quantification: To Be or Not to Be?](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141910)  
*  [paper: The impact of amplification on differential expression analyses by RNA-seq](http://biorxiv.org/content/early/2015/12/28/035493) Computational removal of read duplicates is not recommended for differential expression analysis.
*  [paper: Normalization of RNA-seq data using factor analysis of control genes or samples](https://www.ncbi.nlm.nih.gov/pubmed/25150836 "Risso D et al. Nat Biotechnol 2014"): About spike-ins control and R normalization strategy - remove unwanted variation (RUV).
*  [NVT](https://github.com/Edert/NVT) - an R package for the assessment of RNA-Seq normalization methods.
*  [paper: Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions](http://arxiv.org/abs/1609.00959)  **A must read for RNA-seq normalization and to understand the assumptions for each normalization method!**
* [YARN](http://bioconductor.org/packages/release/bioc/html/yarn.html) bioc package: Robust Multi-Condition RNA-Seq Preprocessing and Normalization.
* [Smooth quantile normalization](https://github.com/stephaniehicks/qsmooth) or qsmooth is a generalization of quantile normalization, which is an average of the two types of assumptions about the data generation process: quantile normalization and quantile normalization between groups.

#### Traditional way of RNA-seq analysis 

* Two nature protocols for RNA-seq analysis  
[Count-based differential expression analysis of RNA sequencing data using R and Bioconductor](http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html)  Based on **DESeq and EdgeR**.  
[Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)  

* A nice tutorial from f1000 research [RNA-Seq workflow: gene-level exploratory analysis and differential expression](http://f1000research.com/articles/4-1070/v1) from Michael Love who is the author of DESeq2.

* [f1000 bioconductor workflow: RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR ](http://f1000research.com/articles/5-1408/v1)
* [f1000 From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline](http://f1000research.com/articles/5-1438/v2) by Gordon Smith.

A post from [Nextgeneseek](http://nextgenseek.com/2015/03/three-papers-on-new-rna-seq-methods-offer-a-new-way-to-do-rna-seq-analysis/)

[QuickRNASeq lifts large-scale RNA-seq data analyses to the next level of automation and interactive visualization](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2356-9)  

>The three papers kind of replaces earlier tools from Salzberg’s group (**Bowtie/TopHat,Cufflinks, and Cuffmerge**)   
they offer a totally new way to go from raw RNA-seq reads to differential expression analysis:  
align RNA-seq reads to genome ([HISAT](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)instead of Bowtie/TopHat, STAR),  
assemble transcripts and estimate expression ([StringTie](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3122.html) instead of Cufflinks), and  
perform differential expression analysis ([Ballgown](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3172.html) instead of Cuffmerge).  

[nature protocol:Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html) 08/11/2016

[Simulation-based comprehensive benchmarking of RNA-seq aligners](http://www.nature.com/nmeth/journal/vaop/ncurrent/pdf/nmeth.4106.pdf)
A nature method paper.  
>We found that performance varied by genome complexity, and accuracy and popularity were poorly correlated. The most widely cited tool underperforms for most metrics, particularly when using default settings

[**RapMap**](http://biorxiv.org/content/early/2015/10/28/029652): A Rapid, Sensitive and Accurate Tool for Mapping RNA-seq Reads to Transcriptomes. From Sailfish group.

* [BitSeq](http://bitseq.github.io/) Transcript isoform level expression and differential expression estimation for RNA-seq
* [Dumpster diving in RNA-sequencing to find the source of every last read](http://biorxiv.org/content/early/2016/05/13/053041) [ROP](https://github.com/smangul1/rop) is a computational protocol aimed to discover the source of all reads, which originated from complex RNA molecules, recombinant antibodies and microbial communities.

**For mapping based methods, usually the raw reads are mapped to transcriptome or genome (need to model gaps by exon-exon junction), and then a gene/transcript level counts are obtained by**:     
* [HTSeq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html): one of the most popular counting tool, but it is slow.
* [featureCounts](http://bioinf.wehi.edu.au/featureCounts/): much faster, use mulitple threads.  
* [VERSE](https://github.com/qinzhu/VERSE): built on `featureCounts`, integrate `HTseq`.  
* [eXpress](http://cdwscience.blogspot.com/2014/02/mrna-quantification-via-express.html). 

Finally, differential expression is carried out by   
* [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [EdgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)  
* [limma Voom](http://www.statsci.org/smyth/pubs/VoomPreprint.pdf) 
* [EBseq](http://www.bioconductor.org/packages/2.14/bioc/html/EBSeq.html) An R package for gene and isoform differential expression analysis of RNA-seq data
* [JunctionSeq](http://hartleys.github.io/JunctionSeq/) differential usage of exons and splice junctions in High-Throughput, Next-Generation RNA-Seq datasets. The methodology is heavily based on the DEXSeq bioconductor package.The core advantage of JunctionSeq over other similar tools is that it provides a powerful automated tools for generating readable and interpretable plots and tables to facilitate the interpretation of the results. An example results report is available [here](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/exampleResults/testForDU.html).

* [MetaSeq](http://bioconductor.org/packages/2.13/bioc/html/metaSeq.html) Meta-analysis of RNA-Seq count data in multiple studies  
* [derfinder](http://www.bioconductor.org/packages/release/bioc/html/derfinder.html) Annotation-agnostic differential expression analysis of RNA-seq data at base-pair resolution  
* [DGEclust](http://dvav.me/dgeclust/) is a program for clustering and differential expression analysis of expression data generated by next-generation sequencing assays, such as RNA-seq, CAGE and others
* [Degust](http://vicbioinformatics.com/degust/index.html): Perform RNA-seq analysis and visualisation. Simply upload a CSV file of read counts for each replicate; then view your DGE data.
* [Vennt](http://drpowell.github.io/vennt/) Dynamic Venn diagrams for Differential Gene Expression.
* [Glimma](http://bioconductor.org/packages/release/bioc/html/Glimma.html)Interactive HTML graphics for RNA-seq data.
#### Extra Notes

* [In RNA-Seq, 2 != 2: Between-sample normalization](https://haroldpimentel.wordpress.com/2014/12/08/in-rna-seq-2-2-between-sample-normalization/)  
* [RPKM/FPKM, TPM and raw counts for RNA-seq](http://crazyhottommy.blogspot.com/2015/06/rpkmfpkm-tpm-and-raw-counts-for-rna-seq.html)  
* [Youtube video counts vs TPM](https://www.youtube.com/watch?v=ztyjiCCt_lM)

### Benchmarking 
[bcbio.rnaseq](https://github.com/roryk/bcbio.rnaseq)    
[RNAseqGUI](http://bioinfo.na.iac.cnr.it/RNASeqGUI/Manual.html). I have used several times. looks good.  
[compcodeR](http://bcf.isb-sib.ch/data/compcodeR/)    
[paper: Benchmark Analysis of Algorithms for Determining and Quantifying Full-length mRNA Splice Forms from RNA-Seq Data](http://bioinformatics.oxfordjournals.org/content/early/2015/09/03/bioinformatics.btv488)  
[paper: Comparative evaluation of isoform-level gene expression estimation algorithms for RNA-seq and exon-array platforms](http://bib.oxfordjournals.org/content/early/2016/02/26/bib.bbw016.short?rss=1)  
[paper:A benchmark for RNA-seq quantification pipelines](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0940-1)

#### Map free 

*  [RNASkim](https://github.com/zzj/RNASkim)
*  [Salmon: Accurate, Versatile and Ultrafast Quantification from RNA-seq Data using Lightweight-Alignment](http://biorxiv.org/content/early/2015/06/27/021592). It is the sucessor of [Salfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/downloads.html)  I have used Salfish once, and it is super-fast! Salmon is supposed to be even better. [tutorial](https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/AGENDA.md)
*  [Kallisto](http://nextgenseek.com/2015/05/kallisto-a-new-ultra-fast-rna-seq-quantitation-method/) from Lior Patcher's lab. [paper: Near-optimal probabilistic RNA-seq quantification](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3519.html) 
*  [sleuth](http://pachterlab.github.io/sleuth/) works with Kallisto for differential expression. 
*  [NORA: A tool for transcript quantification where accuracy matters](https://bioturing.com/nora) claims to be more accurate than RSEM, Salmon and Kallisto.
*  [paper: Differential analysis of RNA-Seq incorporating quantification uncertainty](http://biorxiv.org/content/early/2016/06/10/058164) by Sleuth from Lior Pachter group.
*  [Differential analysis of RNA-Seq incorporating quantification uncertainty: sleuth](http://biorxiv.org/content/early/2016/06/10/058164)
*  [Reanalysis of published RNA-Seq data using kallisto and sleuth](http://lair.berkeley.edu/) based on shiny.
*  [tximport: import and summarize transcript-level estimates for gene-level analysis](https://github.com/mikelove/tximport/blob/master/vignettes/tximport.md) now on [bioconductor](http://bioconductor.org/packages/devel/bioc/html/tximport.html)   
* [tximeta](https://github.com/mikelove/tximeta) by Mike Love. Import transcript quantification into R/Bioconductor with automatic annotation metadata.
*  [f1000 research paper Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v2) from Mike love et.al.
*  [RATS](https://github.com/bartongroup/RATS): Relative Abundance of Transcripts: An R package for the detection of Differential. Transcript isoform Usage.  

>It provides a method to detect changes in the relative abundance of the alternative transcripts (isoforms) of genes. This is called Differential Transcript Usage (DTU).

>Detecting DTU is supplementary to the quantification of transcripts by tools like Salmon, Sailfish and Kallisto and the detection of Differential Transcript Expression (DTE) by tools such as Sleuth.

I particularly like the figure in the tutorial showing the differences among DTU, DTE and DEG.
The paper transcript-level estimates improve gene-level inferences above also talks about the differences:  
1) differential gene expression (DGE) studies, where the overall transcriptional output of each gene is compared between conditions;    2) differential transcript/exon usage (DTU/DEU) studies, where the composition of a gene’s isoform abundance spectrum is compared between conditions, or  
3) differential transcript expression (DTE) studies, where the interest lies in whether individual transcripts show differential expression between conditions.  
![](https://cloud.githubusercontent.com/assets/4106146/18142832/e86c1dd4-6f84-11e6-8efb-8dd404e2942a.png)

* [Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification](https://f1000research.com/articles/7-952/v1)

! [MATS](http://rnaseq-mats.sourceforge.net/rmatsdockerbeta/) is a computational tool to detect differential alternative splicing events from RNA-Seq data. The statistical model of MATS calculates the P-value and false discovery rate that the difference in the isoform ratio of a gene between two conditions exceeds a given user-defined threshold.

### Blog posts on Kallisto/Salmon
1. [Comparing unpublished RNA-Seq gene expression quantifiers](http://nxn.se/post/118321890480/comparing-unpublished-rna-seq-gene-expression)  
2. [Kallisto, a new ultra fast RNA-seq quantitation method](http://nextgenseek.com/2015/05/kallisto-a-new-ultra-fast-rna-seq-quantitation-method/) from Next GEN SEEK
3. [kallisto paper summary: Near-optimal RNA-seq quantification](http://nextgenseek.com/2015/05/kallisto-paper-summary-near-optimal-rna-seq-quantification/) from Next GEN SEEK
4. [Not-quite alignments: Salmon, Kallisto and Efficient Quantification of RNA-Seq data](http://robpatro.com/blog/?p=248)  
5. [Using Kallisto for gene expression analysis of published RNAseq data](https://benchtobioinformatics.wordpress.com/2015/07/10/using-kallisto-for-gene-expression-analysis-of-published-rnaseq-data/)  
6. [How accurate is Kallisto?](http://genomespot.blogspot.com/2015/08/how-accurate-is-kallisto.html) from Mark Ziemann  
7. [ALIGNMENT FREE TRANSCRIPTOME QUANTIFICATION](http://sjcockell.me/2015/05/18/alignment-free-transcriptome-quantification/)
8. [A sleuth for RNA-seq](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)  
9. [Using Salmon, Sailfish and Sleuth for differential expression](http://robpatro.com/blog/?p=291)
10. [Road-testing Kallisto](https://cgatoxford.wordpress.com/2015/10/12/road-testing-kallisto/)
11. [Why you should use alignment-independent quantification for RNA-Seq](https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/)


A biostar [post](https://www.biostars.org/p/143458/#157303): **Do not feed rounded estimates of gene counts from kallisto into DESeq2** (please make sure you read through all the comments, and now there is a suggested workflow for feeding rounded estimates of gene counts to DESeq etc)   
>There is some confusion in the answers to this question that hopefully I can clarify with the three comments below:

>1. kallisto produces estimates of transcript level counts, and therefore to obtain an estimate of the number of reads from a gene the correct thing to do is to sum the estimated counts from the constituent transcripts of that gene. Of note in the language above is the word "estimate", which is necessary because in many cases reads cannot be mapped uniquesly to genes. However insofar as obtaining a good estimate, the approach of kallisto (and before it Cufflinks, RSEM, eXpress and other "transcript level quantification tools") is superior to naïve "counting" approaches for estimating the number of reads originating from a gene. This point has been argued in many papers; among my own papers it is most clearly explained and demonstrated in Trapnell et al.  2013. 

>2. Although estimated counts for a gene can be obtained by summing the estimated counts of the constituent transcripts from tools such as kallisto, and the resulting numbers can be rounded to produce integers that are of the correct format for tools such as DESeq, the numbers produced by such an approach do not satisfy the distributional assumptions made in DESeq and related tools. For example, in DESeq2, counts are modeled "as following a negative binomial distribution". This assumption is not valid when summing estimated counts of transcripts to obtain gene level counts, hence the justified concern of Michael Love that plugging in sums of estimated transcript counts could be problematic for DESeq2. In fact, even the estimated transcript counts themselves are not negative binomial distributed, and therefore also those are not appropriate for plugging into DESeq2. His concern is equally valid with many other "count based" differential expression tools.

>3. Fortunately there is a solution for performing valid statistical testing of differential abundance of individual transcripts, namely the method implemented in sleuth. The approach is described here. To test for differential abundance of genes, one must first address the question of what that means. E.g. is a gene differential if at least one isoform is? or if all the isoforms are? The tests of sleuth are performed at the granularity of transcripts, allowing for downstream analysis that can capture the varied questions that might make biological sense in specific contexts.

>In summary, please do not plug in rounded estimates of gene counts from kallisto into DESeq2 and other tools. While it is technically possible, it is not statistically advisable. Instead, you should use tools that make valid distributional assumptions about the estimates.

**However, Charlotte Soneson, Mike Love and Mark Robinson [showed in a f1000 paper: Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v1) that rounded values from transcript level can be [fed into DESeq2 etc](https://www.biostars.org/p/143458/#144003) for gene-level differential expression, and it is [valid](https://twitter.com/markrobinsonca/status/659823934804766720) and preferable in many ways**.

Thanks [Rob Patro](https://twitter.com/nomad421) for pointing it out!

* [artemis](https://github.com/RamsinghLab/artemis): RNAseq analysis, from raw reads to pathways, typically in a few minutes. Mostly by wrapping `Kallisto` and caching everything we possibly can.
* [isolator](https://github.com/rob-p/isolator):Rapid and robust analysis of RNA-Seq experiments.  

>Isolator has a particular focus on producing stable, consistent estimates. Maximum likelihood approaches produce unstable point estimates: small changes in the data can result in drastically different results, conflating downstream analysis like clustering or PCA. Isolator produces estimates that are in general, simultaneously more stable and more accurate other methods

### Circular RNA

* a bunch of [tools](https://github.com/dieterich-lab/)  from Dieterich Lab in github.

### Batch effects
[TACKLING BATCH EFFECTS AND BIAS IN TRANSCRIPT EXPRESSION](http://mikelove.github.io/eurobioc2015/#/slide-1) by mike love  
[paper:Tackling the widespread and critical impact of batch effects in high-throughput data](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html) by Jeffrey T. Leek in Rafael A. Irizarry's lab.  
[A reanalysis of mouse ENCODE comparative gene expression data](http://f1000research.com/articles/4-121/v1)  
[Is it species or is it batch? They are confounded, so we can't know](http://simplystatistics.org/2015/05/20/is-it-species-or-is-it-batch-they-are-confounded-so-we-cant-know/)    
[Mouse / Human Transcriptomics and Batch Effects](https://rmflight.github.io/posts/2015/06/mouse_human_transcriptomics.html)  
[Meta-analysis of RNA-seq expression data across species, tissues and studies](http://www.genomebiology.com/2015/16/1/287):**Interspecies clustering by tissue is the predominantly observed pattern among various studies under various distance metrics and normalization methods**
[Surrogate Variable Analysis:SVA bioconductor](https://www.bioconductor.org/packages/3.3/bioc/html/sva.html)  
[Paper Summary: Systematic bias and batch effects in single-cell RNA-Seq data](http://nextgenseek.com/2016/01/paper-summary-systematic-bias-and-batch-effects-in-single-cell-rna-seq-data/)  
[Modeling and correcting fragment sequence bias for RNA-seq](https://github.com/mikelove/alpine): alpine bioconductor package from Mike Love.  
[BatchQC](https://github.com/mani2012/BatchQC): interactive software for evaluating sample and batch effects in genomic data.
[A framework for RNA quality correction in differential expression analysis](http://www.biorxiv.org/content/early/2016/09/09/074245)

### Databases
* [BioJupies Automatically Generates RNA-seq Data Analysis Notebooks](https://amp.pharm.mssm.edu/biojupies/) With BioJupies you can produce in seconds a customized, reusable, and interactive report from your own raw or processed RNA-seq data through a simple user interface
* [RNA meta analysis](https://rnama.com/docs/search-evaluation) has ~26,700 studies (5,717 RNA-Seq and 20,955 Microarray). https://rnama.com/  Based on 750 manually labeled studies, our clustering algorithm correctly identifies 91% of sample groups.
* [refine.bio](https://www.refine.bio/)will have harmonized over 60,000 gene expression experiments
* [ReCount is an online resource consisting of RNA-seq gene count datasets built using the raw data from 18 different studies](http://bowtie-bio.sourceforge.net/recount/) updated version [here](https://jhubiostatistics.shinyapps.io/recount/)
* [Recount2-FANTOM](https://www.biorxiv.org/content/10.1101/659490v1) Recounting the FANTOM Cage Associated Transcriptome. Long non-coding RNAs (lncRNAs.
* [The conquer (consistent quantification of external rna-seq data)](http://imlspenticton.uzh.ch:3838/conquer/) repository is developed by Charlotte Soneson and Mark D Robinson at the University of Zurich, Switzerland. single cell RNA-seq data sets.
* [The Lair](http://pachterlab.github.io/lair/): a resource for exploratory analysis of published RNA-Seq data. From Lior Pachter group!
* [The Digital Expression Explorer(dee)](http://dee.bakeridi.edu.au/index.html) The Digital Expression Explorer (DEE) is a repository of digital gene expression profiles mined from public RNA-seq data sets. These data are obtained from NCBI Short Read Archive.  
[blog post for it](http://genomespot.blogspot.com/2015/10/introducing-digital-expression-explorer.html?utm_content=bufferb6214&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)
* [dee2](http://dee2.io/) Digital Expression Explorer 2. Digital Expression Explorer 2 (DEE2) is a repository of uniformly processed RNA-seq data mined from public data obtained from NCBI Short Read Archive. By Ziemann Mark et.al! Version 2 of dee.
* [Extracting allelic read counts from 250,000 human sequencing runs in Sequence Read Archive](https://www.biorxiv.org/content/early/2018/08/07/386441?rss=1) [data](https://www.synapse.org/#!Synapse:syn11415602/wiki/492470)
* [SHARQ Search public, human, RNA-seq experiments by cell, tissue type, and other features | Indexing 19807 files](http://sharq.compbio.cs.cmu.edu/)
* [MetaSRA: normalized sample-specific metadata for the Sequence Read Archive](http://biorxiv.org/content/early/2016/11/30/090506)
* [ARCHS4: Massive Mining of Publicly Available RNA-seq Data from Human and Mouse](https://amp.pharm.mssm.edu/archs4/) ARCHS4 provides access to gene counts from HiSeq 2000, HiSeq 2500 and NextSeq 500 platforms for human and mouse experiments from GEO and SRA. 
* [iDEP-reads: Uniformlly processed public RNA-Seq data](http://bioinformatics.sdstate.edu/reads/) Read counts data for 5,470 human and mouse datasets from ARCHS4 v6 and 12,670 datasets from DEE2 for 9 model organisms by steven Ge.
* [SRA-explorer](https://ewels.github.io/sra-explorer/) This tool aims to make datasets within the Sequence Read Archive more accessible.
* [OmicIDX on BigQuery](https://seandavi.github.io/2019/06/omicidx-on-bigquery/) by Sean Davis who develped SRAdb at NIH. In practice, the OmicIDX mines data from the NCBI Sequence Read Archive (SRA) and NCBI Biosample databases (updated daily).
* [RESTful RNA-seq Analysis API](http://www.ebi.ac.uk/about/news/service-news/new-restful-rna-seq-analysis-api) A simple RESTful API to access analysis results of all public RNAseq data for nearly 200 species in European Nucleotide Archive.
* [intropolis](https://github.com/nellore/intropolis) is a list of exon-exon junctions found across **21,504** human RNA-seq samples on the Sequence Read Archive (SRA) from spliced read alignment to hg19 with Rail-RNA. Two files are provided:
* [ExpressionAtlas bioconductor package](http://www.bioconductor.org/packages/release/bioc/html/ExpressionAtlas.html):
>This package is for searching for datasets in EMBL-EBI Expression Atlas, and downloading them into R for further analysis. Each Expression Atlas dataset is represented as a SimpleList object with one element per platform. Sequencing data is contained in a SummarizedExperiment object, while microarray data is contained in an ExpressionSet or MAList object.
* [GTEx Resources in the UCSC Browser](http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=gtexGene) [signal track on trackhub](http://genome.ucsc.edu/cgi-bin/hgHubConnect)
* [batch recompute ~20,000 RNA-seq samples from larget sequencing project such as TCGA, TARGET and GETEX](https://genome-cancer.soe.ucsc.edu/proj/site/xena/datapages/?host=https://toil.xenahubs.net). Used `hg38` and `gencode v21` as annotation.
* [A cloud-based workflow to quantify transcript-expression levels in public cancer compendia](http://biorxiv.org/content/early/2016/07/12/063552) used kallisto for TCGA/CCLE datasets and gencode v24 as annotation.
* [OMics Compendia Commons (OMiCC)](https://omicc.niaid.nih.gov/) OMiCC is a community-based, biologist-friendly web platform for creating and (meta-) analyzing annotated gene-expression data compendia across studies and technology platforms for more than 24,000 human and mouse studies from Gene Expression Omnibus (GEO)
* [GEIPA](http://gepia.cancer-pku.cn/) interactively explore TCGA expression data, survival etc
* [GEOdiver](http://www.geodiver.co.uk/) An easy to use web tool for analysing GEO datasets.
* [ScanGEO](http://scangeo.dartmouth.edu/ScanGEO/) - parallel mining of high-throughput gene expression data
* [shinyGEO](http://gdancik.github.io/shinyGEO/) a web-based application for performing differential expression and survival analysis on Gene Expression Omnibus datasets.
* [GREIN](https://hub.docker.com/r/ucbd2k/grein/): An interactive web platform for re-analyzing GEO RNA-seq data
* [ImaGEO](http://bioinfo.genyo.es/imageo/) Integrative Meta-Analysis of GEO Data.
* [Expression Atlas update--an integrated database of gene and protein expression in humans, animals and plants](http://www.ebi.ac.uk/gxa/home) It consists of selected microarray and RNA-sequencing studies from ArrayExpress, which have been manually curated, annotated with ontology terms, checked for high quality and processed using standardised analysis methods. Since the last update, Atlas has grown seven-fold (1572 studies as of August 2015), and incorporates baseline expression profiles of tissues from Human Protein Atlas, GTEx and FANTOM5, and of cancer cell lines from ENCODE, CCLE and Genentech projects.
* [DCTD Releases A New Resource for Exploring Cell Line Transcriptional Responses to Anti-Cancer Agents: The NCI Transcriptional Pharmacodynamics Workbench](https://tpwb.nci.nih.gov/GeneExpressionNCI60/TPWorkbench/)
* [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb/) a database for gene expression profiling in human single cell by RNA-seq
* [JingleBells](http://jinglebells.bgu.ac.il/) - A repository of standardized single cell RNA-Seq datasets for analysis and visualization in IGV of the raw reads at the single cell level. Currently focused on immune cells. (http://www.jimmunol.org/content/198/9/3375.long)


### Gene Set enrichment analysis
* [Gene set analysis approaches for RNA-seq data: performance evaluation and application guideline](http://bib.oxfordjournals.org/content/early/2015/09/04/bib.bbv069.long)  
* [Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview](http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html)
* [Metascape](http://metascape.org/gp/index.html) a web server for gene-set analysis.
* [GSEA](http://software.broadinstitute.org/gsea/index.jsp) from Broad Institute.

### Pathway analysis
* [Statistical analysis and visualization of functional profiles for gene and gene clusters: bioconductor clusterProfiler](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html) by GuangChuang Yu from University of HongKong. Can do many jobs and GSEA like figure. It is very useful and I will give it a try besides 
* [GAGE](http://bioconductor.org/packages/release/bioc/html/gage.html).  
* [DAVID](https://david.ncifcrf.gov/):The Database for Annotation, Visualization and Integrated Discovery (DAVID ). *UPDATED in 2016!!!*
* [EGSEA](https://bioconductor.org/packages/release/bioc/html/EGSEA.html)
 Ensemble of Gene Set Enrichment Analyses. By Gordon Smith. take a look!
* [DESeq to fgsea](https://stephenturner.github.io/deseq-to-fgsea/) tutorial by Stephen Turner.
* [Lightweight Iterative Gene set Enrichment in R (LIGER)](https://github.com/JEFworks/liger) by Jean Fan.
* [Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap](https://www.nature.com/articles/s41596-018-0103-9) 2019 Nature Protocol

### Fusion gene detection
* [arriba](https://github.com/suhrig/arriba)Fast and accurate gene fusion detection from RNA-Seq data. top performer of ICGC-TCGA DREAM competition.
* [fusioncatcher](https://github.com/ndaniel/fusioncatcher)  
* [PRADA](https://github.com/crazyhottommy/PRADA_pipeline_Verhaak_lab) from our lab
* [Fusion Matcher](https://github.com/ErasmusMC-Bioinformatics/fuma): Match predicted fusions according to chromosomal location or gene annotation(s)
* [paper:Comprehensive evaluation of fusion transcript detection algorithms and a meta-caller to combine top performing methods in paired-end RNA-seq data](https://nar.oxfordjournals.org/content/early/2015/11/17/nar.gkv1234.full)  
* [paper: Comparative assessment of methods for the fusion transcripts detection from RNA-Seq data](http://www.nature.com/articles/srep21597)
* [chimera](https://www.bioconductor.org/packages/release/bioc/html/chimera.html) A package for secondary analysis of fusion products.
* [Pegasus](https://sourceforge.net/p/pegasus-fus/wiki/Main%20Manual/) Fusion Annotation and Prediction.
* [Oncofuse](https://github.com/mikessh/oncofuse) is a framework designed to estimate the oncogenic potential of de-novo discovered gene fusions. It uses several hallmark features and employs a bayesian classifier to provide the probability of a given gene fusion being a driver mutation.
* [chimeraviz](http://bioconductor.org/packages/devel/bioc/html/chimeraviz.html) Visualization tools for gene fusions.
* [SQUID](https://github.com/Kingsford-Group/squid): [Transcriptomic Structural Variation Detection from RNA-seq](https://www.biorxiv.org/content/early/2017/09/06/162776)

### Alternative splicing
* [SplicePlot: a tool for visualizing alternative splicing](http://montgomerylab.stanford.edu/spliceplot/index.html) Sashimi plots
* [Multivariate Analysis of Transcript Splicing (MATS)](http://rnaseq-mats.sourceforge.net/)
* [SNPlice](https://code.google.com/p/snplice/) is a software tool to find and evaluate the co-occurrence of single-nucleotide-polymorphisms (SNP) and altered splicing in next-gen mRNA sequence reads. SNPlice requires, as input: genome aligned reads, exon-intron-exon junctions, and SNPs. exon-intron-exon junctions and SNPs may be derived from the reads directly, using, for example, TopHat2 and samtools, or they may be derived from independent sources
* [Visualizing Alternative Splicing](http://vials.io/vials/) [github page](https://github.com/Caleydo/vials)
* [spladder](https://github.com/ratschlab/spladder) Tool for the detection and quantification of alternative splicing events from RNA-Seq data
* [SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa) This tool generates different Alternative Splicing (AS) events and calculates the PSI ("Percentage Spliced In") value for each event exploiting the fast quantification of transcript abundances from multiple samples.
* [IntSplice](http://www.med.nagoya-u.ac.jp/neurogenetics/IntSplice/): Upload a VCF (variant call format) file to predict if an SNV (single nucleotide variation) from intronic positions -50 to -3 is pathogenic or not.
* [Whippet](http://www.biorxiv.org/content/early/2017/07/03/158519): an efficient method for the detection and quantification of alternative splicing reveals extensive transcriptomic complexity


### microRNAs and non-coding RNAs
* [miARma-Seq workflow](http://miarmaseq.cbbio.es/) miRNA-Seq And RNA-Seq Multiprocess Analysis tool, a comprehensive pipeline analysis suite designed for mRNA, miRNA and circRNA identification and differential expression analysis, applicable to any sequenced organism.
* [All the tools you need to analyse your miRNAs:[tools4miRNAs](http://tools4mirs.org/)
* paper [Evaluation of microRNA alignment techniques](http://rnajournal.cshlp.org/content/early/2016/06/09/rna.055509.115?top=1)
* [protocol: Analysis RNA-seq and Noncoding RNA](http://link.springer.com/protocol/10.1007%2F978-1-4939-6380-5_11)

### transcriptional pausing 

* GRO-seq
* RNApol2 ChIP-seq
* [iRNA-seq: computational method for genome-wide assessment of acute transcriptional regulation from total RNA-seq data](http://nar.oxfordjournals.org/content/43/6/e40.long)  

### intron retention
* [IRFinder](https://github.com/williamritchie/IRFinder)

### Allel specific expression
* paper [Tools and best practices for data processing in allelic expression analysis](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
* paper [Characterizing noise structure in single-cell RNA-seq distinguishes genuine from technical stochastic allelic expression](http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html)  
* paper [Tools and best practices for data processing in allelic expression analysis](http://www.genomebiology.com/2015/16/1/195)

### immnune related
* [ImReP](https://sergheimangul.wordpress.com/imrep/) is a computational method for rapid and accurate profiling of the adaptive immune repertoire from regular RNA-Seq data.
* [Comprehensive analyses of tumor immunity: implications for cancer immunotherapy](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1028-7) 
* [pVAC-Seq](https://github.com/griffithlab/pVAC-Seq) is a cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVAC-Seq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. 
* [JingleBells] (http://jinglebells.bgu.ac.il/) - A repository of standardized single cell RNA-Seq datasets for analysis and visualization in IGV of the raw reads at the single cell level. Currently focused on immune cells. (http://www.jimmunol.org/content/198/9/3375.long)
* [immunedeconv](https://grst.github.io/immunedeconv/) - an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data.

### Reads from xenografts

* [Xenosplit](https://github.com/goknurginer/XenoSplit) XenoSplit is a fast computational tool to detect the true origin of the graft RNA-Seq and DNA-Seq libraries prior to profiling of patient-derived xenografts (PDXs). 

### single cell tutorials
* [Course material in notebook format for learning about single cell bioinformatics methods](https://github.com/YeoLab/single-cell-bioinformatics)
* [Analysis of single cell RNA-seq data course, Cambridge University](https://github.com/hemberg-lab/scRNA.seq.course) Great tutorial!
* [f1000 workflow paper A step-by-step workflow for low-level analysis of single-cell RNA-seq data](http://f1000research.com/articles/5-2122/v1) by Aaron Lun, the athour of diffHiC, GenomicInteractions and csaw.
* [2016 Bioconductor workshop: Analysis of single-cell RNA-seq data with R and Bioconductor](https://github.com/drisso/bioc2016singlecell)
* [paper: Single-Cell Transcriptomics Bioinformatics and Computational Challenges](http://journal.frontiersin.org/article/10.3389/fgene.2016.00163/full)


### single cell RNA-seq normalization
* [paper: Assessment of single cell RNA-seq normalization methods](http://biorxiv.org/content/early/2016/07/17/064329)
* [paper: A practical guide to single-cell RNA-sequencing for biomedical research and clinical applications](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0467-4)
* [Normalizing single-cell RNA sequencing data: challenges and opportunities](https://www.nature.com/nmeth/journal/v14/n6/full/nmeth.4292.html) Nature Methods
* [SinQC: A Method and Tool to Control Single-cell RNA-seq Data Quality](http://www.morgridge.net/SinQC.html).
* [Scone](https://github.com/YosefLab/scone) Single-Cell Overview of Normalized Expression data

### single cell batch effect
* [Overcoming confounding plate effects in differential expression analyses of single-cell RNA-seq data](http://biorxiv.org/content/early/2016/09/08/073973)

### Single cell RNA-seq

* [a collection of single RNA-seq tools by Sean Davis ](https://github.com/seandavi/awesome-single-cell)
* [paper: Design and computational analysis of single-cell RNA-sequencing experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)
* [paper by Mark Robinson: Bias, Robustness And Scalability In Differential Expression Analysis Of Single-Cell RNA-Seq Data](http://biorxiv.org/content/early/2017/05/28/143289)

> Considerable differences are found between the methods in terms of the number and characteristics of the genes that are called differentially expressed. Pre-filtering of lowly expressed genes can have important effects on the results, particularly for some of the methods originally developed for analysis of bulk RNA-seq data. Generally, however, **methods developed for bulk RNA-seq analysis do not perform notably worse than those developed specifically for scRNA-seq.**

* [paper: Power Analysis of Single Cell RNA‐Sequencing Experiments](http://biorxiv.org/content/early/2016/09/08/073692)
* [paper: The contribution of cell cycle to heterogeneity in single-cell RNA-seq data](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3498.html)
* [paper: Batch effects and the effective design of single-cell gene expression studies](http://biorxiv.org/content/early/2016/07/08/062919)
* [On the widespread and critical impact of systematic bias and batch effects in single-cell RNA-Seq data](http://biorxiv.org/content/early/2015/08/25/025528) 
* [paper: Comparison of methods to detect differentially expressed genes between single-cell populations](http://www.ncbi.nlm.nih.gov/pubmed/27373736)
* [review: Single-cell genome sequencing: current state of the science](http://www.nature.com/nrg/journal/vaop/ncurrent/full/nrg.2015.16.html)  
* [Ginkgo](http://qb.cshl.edu/ginkgo/?q=/ESjKTTeZIdnoGwEB4WTu) A web tool for analyzing single-cell sequencing data.
* [SingleCellExperiment bioc package](http://bioconductor.org/packages/devel/bioc/html/SingleCellExperiment.html) Defines a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.
* [ASAP](http://biorxiv.org/content/early/2016/12/22/096222): a Web-based platform for the analysis and inter-active visualization of single-cell RNA-seq data
* [FASTGenomics](https://fastgenomics.org), an online platform to share single-cell RNA sequencing data and perform analyses using reproducible workflows. Users can upload their own data and use standard or customized workflows for the exploration and analysis of gene expression data ([Scholz et al. 2018](https://doi.org/10.1101/272476)).
* [Seurat](http://www.satijalab.org/seurat.html) is an R package designed for the analysis and visualization of single cell RNA-seq data. It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
* [R package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data](http://bioconductor.org/packages/devel/bioc/html/sincell.html)  
* [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) Differential expression and time-series analysis for single-cell RNA-Seq and qPCR experiments.
* Single Cell Differential Expression: bioconductor package [scde](http://bioconductor.org/packages/devel/bioc/html/scde.html)
* [Sincera](https://research.cchmc.org/pbge/sincera.html):A Computational Pipeline for Single Cell RNA-Seq Profiling Analysis. Bioconductor package will be available soon. 
* [MAST](https://github.com/RGLab/MAST): a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data
* [scDD](http://biorxiv.org/content/early/2015/11/13/031021): A statistical approach for identifying differential distributions in single-cell RNA-seq experiments
* [Inference and visualisation of Single-Cell RNA-seq Data data as a hierarchical tree structure: bioconductor CellTree](http://bioconductor.org/packages/devel/bioc/html/cellTree.html)  
* [Fast and accurate single-cell RNA-Seq analysis by clustering of transcript-compatibility counts](http://biorxiv.org/content/early/2016/01/15/036863) by Lior Pachter et.al
* [cellity](https://github.com/ti243/cellity): Classification of low quality cells in scRNA-seq data using R.
* [bioconductor: using scran to perform basic analyses of single-cell RNA-seq data](http://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html)
* [scater](https://github.com/davismcc/scater): single-cell analysis toolkit for expression with R
* [Monovar](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3835.html): single-nucleotide variant detection in single cells
* [paper: Comparison of methods to detect differentially expressed genes between single-cell populations](http://m.bib.oxfordjournals.org/content/early/2016/07/02/bib.bbw057.full)
* [Single-cell mRNA quantification and differential analysis with Census](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4150.html)
* [CIDR](https://github.com/VCCRI/CIDR): Ultrafast and accurate clustering through imputation for single-cell RNA-seq data
* [CellView](http://biorxiv.org/content/early/2017/04/04/123810): Interactive Exploration Of High Dimensional Single Cell RNA-Seq Data

### single cell RNA-seq clustering

* [Geometry of the Gene Expression Space of Individual Cells](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004224)
* [pcaReduce](http://biorxiv.org/content/early/2015/09/08/026385): Hierarchical Clustering of Single Cell Transcriptional Profiles.
* [CountClust](https://www.bioconductor.org/packages/3.3/bioc/html/CountClust.html): Clustering and Visualizing RNA-Seq Expression Data using Grade of Membership Models. Fits grade of membership models (GoM, also known as admixture models) to cluster RNA-seq gene expression count data, identifies characteristic genes driving cluster memberships, and provides a visual summary of the cluster memberships
* [FastProject](http://biorxiv.org/content/early/2016/03/12/043463): A Tool for Low-Dimensional Analysis of Single-Cell RNA-Seq Data
* [SNN-Cliq](http://bioinformatics.oxfordjournals.org/content/31/12/1974.full) Identification of cell types from single-cell transcriptomes using a novel clustering method
* [Compare clusterings for single-cell sequencing](http://bioconductor.org/packages/devel/bioc/html/clusterExperiment.html) bioconductor package.The goal of this package is to encourage the user to try many different clustering algorithms in one package structure. We give tools for running many different clusterings and choices of parameters. We also provide visualization to compare many different clusterings and algorithm tools to find common shared clustering patterns.
* [CIDR: Ultrafast and accurate clustering through imputation for single cell RNA-Seq data](http://biorxiv.org/content/early/2016/08/10/068775)
* [SC3](http://bioconductor.org/packages/release/bioc/html/SC3.html)- consensus clustering of single-cell RNA-Seq data.  SC3 achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. Tests on twelve published datasets show that SC3 outperforms five existing methods while remaining scalable, as shown by the analysis of a large dataset containing 44,808 cells. Moreover, an interactive graphical implementation makes SC3 accessible to a wide audience of users, and SC3 aids biological interpretation by identifying marker genes, differentially expressed genes and outlier cells.

### advance of scRNA-seq tech
* [Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding](http://science.sciencemag.org/content/360/6385/176) no isolation of single cells needed!
