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

### RNA-seq specific 

*  [Introduction to RNA-seq analysis youtube video](https://www.youtube.com/watch?v=OEbjHPk20C0&feature=youtu.be&a)  
*  [RNAseq differential expression analysis – NGS2015](https://monsterbashseq.wordpress.com/2015/08/26/rnaseq-differential-expression-analysis-ngs2015/)  
*  [Kallisto and sleuth tutorial](http://pachterlab.github.io/sleuth/starting.html) blazing fast RNA-seq analysis by Lior Patcher's lab.    [A sleuth for RNA-Seq](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)  
*  pathway analysis using [GAGE](https://github.com/ajwije/150826_pathway_analysis/blob/master/Tutorial_150827.Rmd) 
*  [Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview](http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html)
*  [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/14/6/671.full)
*  [RNA-seq tutorial wiki](https://github.com/crazyhottommy/rnaseq_tutorial) Informatics for RNA-seq: A web resource for analysis on the cloud.  
*  [RNA-seqlopedia](http://rnaseq.uoregon.edu/)  Great introduction of RNA-seq from sample preparation to data analysis
*  [RNAseq data analysis from data carpentry](https://github.com/datacarpentry/rnaseq-data-analysis)
*  [paper: Isoform prefiltering improves performance of count-based methods for analysis of differential transcript usage](http://www.genomebiology.com/2016/17/1/12)
*  [paper: A survey of best practices for RNA-seq data analysis](http://www.genomebiology.com/2016/17/1/13)
*  [paper: Cross-platform normalization of microarray and RNA-seq data for machine learning applications](https://peerj.com/articles/1621/#results|discussion). [Tool](https://github.com/greenelab/TDM) 
*  [review: Translating RNA sequencing into clinical diagnostics: opportunities and challenges](http://www.nature.com/nrg/journal/v17/n5/full/nrg.2016.10.html)

### RNA-seq experimental design 
* [Tutorial: Rna Seq Experimental Design For Measuring Differential Gene Expression from biostars](https://www.biostars.org/p/65824/)
* [Scotty - Power Analysis for RNA Seq Experiments](http://bioinformatics.bc.edu/marthlab/scotty/scotty.php)  
* [Experimental Design in Differential Abundance analysis web server](http://edda.gis.a-star.edu.sg/)  
* [Experimental Design in Differential Abundance analysis bioconductor package](http://www.bioconductor.org/packages/devel/bioc/html/EDDA.html)

### Quality Control

* [QoRTs](http://hartleys.github.io/QoRTs/): a comprehensive toolset for quality control and data processing of RNA-Seq experiments  
* [QUaCRS](http://bioserv.mps.ohio-state.edu/QuaCRS/index.php/pages/view/downloads)    
* [RSeQC](http://rseqc.sourceforge.net/) RNA-seq data QC  

### Normalization, quantification, and differential expression

*  [A Comparison of Methods: Normalizing High-Throughput RNA Sequencing Data](http://biorxiv.org/content/early/2015/09/03/026062)
*  [Errors in RNA-Seq quantification affect genes of relevance to human disease](http://www.genomebiology.com/2015/16/1/177)  
*  [A comprehensive evaluation of ensembl, RefSeq, and UCSC annotations in the context of RNA-seq read mapping and gene quantification](http://www.biomedcentral.com/1471-2164/16/97)  
*  [Comparing the normalization methods for the differential analysis of Illumina high-throughput RNA-Seq data](http://www.biomedcentral.com/1471-2105/16/347)
*  paper: [Union Exon Based Approach for RNA-Seq Gene Quantification: To Be or Not to Be?](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141910)  
*  [paper: The impact of amplification on differential expression analyses by RNA-seq](http://biorxiv.org/content/early/2015/12/28/035493) Computational removal of read duplicates is not recommended for differential expression analysis.

#### Traditional way of RNA-seq analysis 

* Two nature protocols for RNA-seq analysis  
[Count-based differential expression analysis of RNA sequencing data using R and Bioconductor](http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html)  Based on **DESeq and EdgeR**.  
[Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)  

A nice tutorial from f1000 research [RNA-Seq workflow: gene-level exploratory analysis and differential expression](http://f1000research.com/articles/4-1070/v1) from Michael Love who is the author of DESeq2.

A post from [Nextgeneseek](http://nextgenseek.com/2015/03/three-papers-on-new-rna-seq-methods-offer-a-new-way-to-do-rna-seq-analysis/)

[QuickRNASeq lifts large-scale RNA-seq data analyses to the next level of automation and interactive visualization](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2356-9)  

>The three papers kind of replaces earlier tools from Salzberg’s group (**Bowtie/TopHat,Cufflinks, and Cuffmerge**)   
they offer a totally new way to go from raw RNA-seq reads to differential expression analysis:  
align RNA-seq reads to genome ([HISAT](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)instead of Bowtie/TopHat, STAR),  
assemble transcripts and estimate expression ([StringTie](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3122.html) instead of Cufflinks), and  
perform differential expression analysis ([Ballgown](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3172.html) instead of Cuffmerge).  

[**RapMap**](http://biorxiv.org/content/early/2015/10/28/029652): A Rapid, Sensitive and Accurate Tool for Mapping RNA-seq Reads to Transcriptomes. From Sailfish group.

* [BitSeq](http://bitseq.github.io/) Transcript isoform level expression and differential expression estimation for RNA-seq

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
*  [Reanalysis of published RNA-Seq data using kallisto and sleuth](http://lair.berkeley.edu/) based on shiny.
*  [tximport: import and summarize transcript-level estimates for gene-level analysis](https://github.com/mikelove/tximport/blob/master/vignettes/tximport.md) now on [bioconductor](http://bioconductor.org/packages/devel/bioc/html/tximport.html)   
*  [f1000 research paper Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v1) from Mike love et.al.


### Blog posts on Kallisto
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


A biostar [post](https://www.biostars.org/p/143458/#157303): **Do not feed rounded estimates of gene counts from kallisto into DESeq2** (please make sure you read through all the comments, and now there is a suggested workflow for feeding rounded estimates of gene counts to DESeq etc)   
>There is some confusion in the answers to this question that hopefully I can clarify with the three comments below:

>1. kallisto produces estimates of transcript level counts, and therefore to obtain an estimate of the number of reads from a gene the correct thing to do is to sum the estimated counts from the constituent transcripts of that gene. Of note in the language above is the word "estimate", which is necessary because in many cases reads cannot be mapped uniquely to genes. However insofar as obtaining a good estimate, the approach of kallisto (and before it Cufflinks, RSEM, eXpress and other "transcript level quantification tools") is superior to naïve "counting" approaches for estimating the number of reads originating from a gene. This point has been argued in many papers; among my own papers it is most clearly explained and demonstrated in Trapnell et al.  2013. 

>2. Although estimated counts for a gene can be obtained by summing the estimated counts of the constituent transcripts from tools such as kallisto, and the resulting numbers can be rounded to produce integers that are of the correct format for tools such as DESeq, the numbers produced by such an approach do not satisfy the distributional assumptions made in DESeq and related tools. For example, in DESeq2, counts are modeled "as following a negative binomial distribution". This assumption is not valid when summing estimated counts of transcripts to obtain gene level counts, hence the justified concern of Michael Love that plugging in sums of estimated transcript counts could be problematic for DESeq2. In fact, even the estimated transcript counts themselves are not negative binomial distributed, and therefore also those are not appropriate for plugging into DESeq2. His concern is equally valid with many other "count based" differential expression tools.

>3. Fortunately there is a solution for performing valid statistical testing of differential abundance of individual transcripts, namely the method implemented in sleuth. The approach is described here. To test for differential abundance of genes, one must first address the question of what that means. E.g. is a gene differential if at least one isoform is? or if all the isoforms are? The tests of sleuth are performed at the granularity of transcripts, allowing for downstream analysis that can capture the varied questions that might make biological sense in specific contexts.

>In summary, please do not plug in rounded estimates of gene counts from kallisto into DESeq2 and other tools. While it is technically possible, it is not statistically advisable. Instead, you should use tools that make valid distributional assumptions about the estimates.

**However, Charlotte Soneson, Mike Love and Mark Robinson [showed in a f1000 paper: Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v1) that rounded values from transcript level can be [fed into DESeq2 etc](https://www.biostars.org/p/143458/#144003) for gene-level differential expression, and it is [valid](https://twitter.com/markrobinsonca/status/659823934804766720) and preferable in many ways**.

Thanks [Rob Patro](https://twitter.com/nomad421) for pointing it out!

* [artemis](https://github.com/RamsinghLab/artemis): RNAseq analysis, from raw reads to pathways, typically in a few minutes. Mostly by wrapping `Kallisto` and caching everything we possibly can.
* [isolator](https://github.com/rob-p/isolator):Rapid and robust analysis of RNA-Seq experiments.  

>Isolator has a particular focus on producing stable, consistent estimates. Maximum likelihood approaches produce unstable point estimates: small changes in the data can result in drastically different results, conflating downstream analysis like clustering or PCA. Isolator produces estimates that are in general, simultaneously more stable and more accurate other methods

### Batch effects
[TACKLING BATCH EFFECTS AND BIAS IN TRANSCRIPT EXPRESSION](http://mikelove.github.io/eurobioc2015/#/slide-1) by mike love  
[paper:Tackling the widespread and critical impact of batch effects in high-throughput data](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html) by Jeffrey T. Leek in Rafael A. Irizarry's lab.  
[A reanalysis of mouse ENCODE comparative gene expression data](http://f1000research.com/articles/4-121/v1)  
[Is it species or is it batch? They are confounded, so we can't know](http://simplystatistics.org/2015/05/20/is-it-species-or-is-it-batch-they-are-confounded-so-we-cant-know/)    
[Mouse / Human Transcriptomics and Batch Effects](https://rmflight.github.io/posts/2015/06/mouse_human_transcriptomics.html)  
[Meta-analysis of RNA-seq expression data across species, tissues and studies](http://www.genomebiology.com/2015/16/1/287):**Interspecies clustering by tissue is the predominantly observed pattern among various studies under various distance metrics and normalization methods**
[Surrogate Variable Analysis:SVA bioconductor](https://www.bioconductor.org/packages/3.3/bioc/html/sva.html)  
[Paper Summary: Systematic bias and batch effects in single-cell RNA-Seq data](http://nextgenseek.com/2016/01/paper-summary-systematic-bias-and-batch-effects-in-single-cell-rna-seq-data/)  

### Databases
* [ReCount is an online resource consisting of RNA-seq gene count datasets built using the raw data from 18 different studies](http://bowtie-bio.sourceforge.net/recount/)
* [The Digital Expression Explorer](http://dee.bakeridi.edu.au/index.html) The Digital Expression Explorer (DEE) is a repository of digital gene expression profiles mined from public RNA-seq data sets. These data are obtained from NCBI Short Read Archive.  
[blog post for it](http://genomespot.blogspot.com/2015/10/introducing-digital-expression-explorer.html?utm_content=bufferb6214&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)
* [SHARQ Search public, human, RNA-seq experiments by cell, tissue type, and other features | Indexing 19807 files](http://sharq.compbio.cs.cmu.edu/)
* [RESTful RNA-seq Analysis API](http://www.ebi.ac.uk/about/news/service-news/new-restful-rna-seq-analysis-api) A simple RESTful API to access analysis results of all public RNAseq data for nearly 200 species in European Nucleotide Archive.
* [intropolis](https://github.com/nellore/intropolis) is a list of exon-exon junctions found across **21,504** human RNA-seq samples on the Sequence Read Archive (SRA) from spliced read alignment to hg19 with Rail-RNA. Two files are provided:
* [ExpressionAtlas bioconductor package](http://www.bioconductor.org/packages/release/bioc/html/ExpressionAtlas.html):
>This package is for searching for datasets in EMBL-EBI Expression Atlas, and downloading them into R for further analysis. Each Expression Atlas dataset is represented as a SimpleList object with one element per platform. Sequencing data is contained in a SummarizedExperiment object, while microarray data is contained in an ExpressionSet or MAList object.

### Gene Set enrichment analysis
* [Gene set analysis approaches for RNA-seq data: performance evaluation and application guideline](http://bib.oxfordjournals.org/content/early/2015/09/04/bib.bbv069.long)  
* [Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview](http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html)  

### Pathway analysis
[Statistical analysis and visualization of functional profiles for gene and gene clusters: bioconductor clusterProfiler](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html) by GuangChuang Yu from University of HongKong. Can do many jobs and GSEA like figure. It is very useful and I will give it a try besides [GAGE](http://bioconductor.org/packages/release/bioc/html/gage.html).
[DAVID](https://david.ncifcrf.gov/):The Database for Annotation, Visualization and Integrated Discovery (DAVID ). *UPDATED in 2016!!!*


### Fusion gene detection
* [fusioncatcher](https://github.com/ndaniel/fusioncatcher)  
* [PRADA](https://github.com/crazyhottommy/PRADA_pipeline_Verhaak_lab) from our lab
* [Fusion Matcher](https://github.com/ErasmusMC-Bioinformatics/fuma): Match predicted fusions according to chromosomal location or gene annotation(s)
* [paper:Comprehensive evaluation of fusion transcript detection algorithms and a meta-caller to combine top performing methods in paired-end RNA-seq data](https://nar.oxfordjournals.org/content/early/2015/11/17/nar.gkv1234.full)  
* [paper: Comparative assessment of methods for the fusion transcripts detection from RNA-Seq data](http://www.nature.com/articles/srep21597)
* [chimera](https://www.bioconductor.org/packages/release/bioc/html/chimera.html) A package for secondary analysis of fusion products.
* [Pegasus](https://sourceforge.net/p/pegasus-fus/wiki/Main%20Manual/) Fusion Annotation and Prediction.
* [Oncofuse](https://github.com/mikessh/oncofuse) is a framework designed to estimate the oncogenic potential of de-novo discovered gene fusions. It uses several hallmark features and employs a bayesian classifier to provide the probability of a given gene fusion being a driver mutation.

### Alternative splicing
* [SplicePlot: a tool for visualizing alternative splicing](http://montgomerylab.stanford.edu/spliceplot/index.html) Sashimi plots
* [Multivariate Analysis of Transcript Splicing (MATS)](http://rnaseq-mats.sourceforge.net/)
* [SNPlice](https://code.google.com/p/snplice/) is a software tool to find and evaluate the co-occurrence of single-nucleotide-polymorphisms (SNP) and altered splicing in next-gen mRNA sequence reads. SNPlice requires, as input: genome aligned reads, exon-intron-exon junctions, and SNPs. exon-intron-exon junctions and SNPs may be derived from the reads directly, using, for example, TopHat2 and samtools, or they may be derived from independent sources
* [Visualizing Alternative Splicing](http://vials.io/vials/) [github page](https://github.com/Caleydo/vials)
* [spladder](https://github.com/ratschlab/spladder) Tool for the detection and quantification of alternative splicing events from RNA-Seq data
* [SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa) This tool generates different Alternative Splicing (AS) events and calculates the PSI ("Percentage Spliced In") value for each event exploiting the fast quantification of transcript abundances from multiple samples

### microRNAs and non-coding RNAs
* [miARma-Seq workflow](http://miarmaseq.cbbio.es/) miRNA-Seq And RNA-Seq Multiprocess Analysis tool, a comprehensive pipeline analysis suite designed for mRNA, miRNA and circRNA identification and differential expression analysis, applicable to any sequenced organism.

### transcriptional pausing 

* GRO-seq
* RNApol2 ChIP-seq
* [iRNA-seq: computational method for genome-wide assessment of acute transcriptional regulation from total RNA-seq data](http://nar.oxfordjournals.org/content/43/6/e40.long)  

### Allel specific expression
* paper [Tools and best practices for data processing in allelic expression analysis](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
* paper [Characterizing noise structure in single-cell RNA-seq distinguishes genuine from technical stochastic allelic expression](http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html)  
* paper [Tools and best practices for data processing in allelic expression analysis](http://www.genomebiology.com/2015/16/1/195)

### Single cell RNA-seq
* [On the widespread and critical impact of systematic bias and batch effects in single-cell RNA-Seq data](http://biorxiv.org/content/early/2015/08/25/025528) 
* [review: Single-cell genome sequencing: current state of the science](http://www.nature.com/nrg/journal/vaop/ncurrent/full/nrg.2015.16.html)  
* [Ginkgo](http://qb.cshl.edu/ginkgo/?q=/ESjKTTeZIdnoGwEB4WTu) A web tool for analyzing single-cell sequencing data.
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
* [Monovar](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3835.html): single-nucleotide variant detection in single cells
