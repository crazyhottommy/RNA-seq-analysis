# RNA-seq analysis
I will put some RNA-seq resources here.  

### General sequencing data analysis materials
* [Next-Gen Sequence Analysis Workshop (2015)](http://angus.readthedocs.org/en/2015/) held by [Titus Brown](http://genomecenter.ucdavis.edu/people/faculty/name/c-titus-brown/)  (now in UC Davis)
* [Fall 2015, BMMB 852: Applied Bioinformatics](http://www.personal.psu.edu/iua1/2015_fall_852/main_2015_fall_852.html) by        [Istvan Albert](http://www.personal.psu.edu/iua1/) from Penn state University. He developed the all-time popular 
   [biostars](https://www.biostars.org/)  
* Steven Turner in UVA is maitaining a list of training opportunities for [genomic data analysis](http://stephenturner.us/edu.html)
*  Jeff Leek group's recommended [genomic papers](https://github.com/jtleek/genomicspapers/)
* [awesome tutorial for NGS file format](http://binf.snipcademy.com/lessons/sequence-file-formats)  

### RNA-seq specific 

*  [Introduction to RNA-seq analysis youtube video](https://www.youtube.com/watch?v=OEbjHPk20C0&feature=youtu.be&a)  
*  [RNAseq differential expression analysis – NGS2015](https://monsterbashseq.wordpress.com/2015/08/26/rnaseq-differential-expression-analysis-ngs2015/)  
*  [Kallisto and sleuth tutorial](http://pachterlab.github.io/sleuth/starting.html) blazing fast RNA-seq analysis by Lior Patcher's lab.    [A sleuth for RNA-Seq](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)  
*  pathway analysis using [GAGE](https://github.com/ajwije/150826_pathway_analysis/blob/master/Tutorial_150827.Rmd)  
*  [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/14/6/671.full)
*  [RNA-seq tutorial wiki](https://github.com/crazyhottommy/rnaseq_tutorial) Informatics for RNA-seq: A web resource for analysis on the cloud.  
*  [RNA-seqlopedia](http://rnaseq.uoregon.edu/)  Great introduction of RNA-seq from sample preparation to data analysis
*  [RNAseq data analysis from data carpentry](https://github.com/datacarpentry/rnaseq-data-analysis)

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

#### Traditional way of RNA-seq analysis 

* Two nature protocols for RNA-seq analysis  
[Count-based differential expression analysis of RNA sequencing data using R and Bioconductor](http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html)  Based on **DESeq and EdgeR**.  
[Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)  

A nice tutorial from f1000 research [RNA-Seq workflow: gene-level exploratory analysis and differential expression](http://f1000research.com/articles/4-1070/v1) from Michael Love who is the author of DESeq2.

A post from [Nextgeneseek](http://nextgenseek.com/2015/03/three-papers-on-new-rna-seq-methods-offer-a-new-way-to-do-rna-seq-analysis/)  

>The three papers kind of replaces earlier tools from Salzberg’s group (**Bowtie/TopHat,Cufflinks, and Cuffmerge**)   
they offer a totally new way to go from raw RNA-seq reads to differential expression analysis:  
align RNA-seq reads to genome ([HISAT](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)instead of Bowtie/TopHat, STAR),  
assemble transcripts and estimate expression ([StringTie](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3122.html) instead of Cufflinks), and  
perform differential expression analysis ([Ballgown](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3172.html) instead of Cuffmerge).  

* [BitSeq](http://bitseq.github.io/) Transcript isoform level expression and differential expression estimation for RNA-seq

**For mapping based methods, usually the raw reads are mapped to transcriptome or genome (need to model gaps by exon-exon junction), and then a gene/transcript level counts are obtained by   
* [HTSeq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)    
* [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)   
* [eXpress](http://cdwscience.blogspot.com/2014/02/mrna-quantification-via-express.html). 

Finally, differential expression is carried out by   
* [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [EdgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)  
* [limma Voom](http://www.statsci.org/smyth/pubs/VoomPreprint.pdf) 
* [EBseq](http://www.bioconductor.org/packages/2.14/bioc/html/EBSeq.html) An R package for gene and isoform differential expression analysis of RNA-seq data

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

#### Map free 

*  [RNASkim](https://github.com/zzj/RNASkim)
*  [Salmon: Accurate, Versatile and Ultrafast Quantification from RNA-seq Data using Lightweight-Alignment](http://biorxiv.org/content/early/2015/06/27/021592). It is the sucessor of [Salfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/downloads.html)  I have used Salfish once, and it is super-fast! Salmon is supposed to be even better.   
*  [Kallisto](http://nextgenseek.com/2015/05/kallisto-a-new-ultra-fast-rna-seq-quantitation-method/) from Lior Patcher's lab.
*  [sleuth](http://pachterlab.github.io/sleuth/) works with Kallisto for differential expression.


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


A biostar [post](https://www.biostars.org/p/143458/#157303): **Do not feed rounded estimates of gene counts from kallisto into DESeq2**   
>There is some confusion in the answers to this question that hopefully I can clarify with the three comments below:

>1. kallisto produces estimates of transcript level counts, and therefore to obtain an estimate of the number of reads from a gene the correct thing to do is to sum the estimated counts from the constituent transcripts of that gene. Of note in the language above is the word "estimate", which is necessary because in many cases reads cannot be mapped uniquely to genes. However insofar as obtaining a good estimate, the approach of kallisto (and before it Cufflinks, RSEM, eXpress and other "transcript level quantification tools") is superior to naïve "counting" approaches for estimating the number of reads originating from a gene. This point has been argued in many papers; among my own papers it is most clearly explained and demonstrated in Trapnell et al.  2013. 

>2. Although estimated counts for a gene can be obtained by summing the estimated counts of the constituent transcripts from tools such as kallisto, and the resulting numbers can be rounded to produce integers that are of the correct format for tools such as DESeq, the numbers produced by such an approach do not satisfy the distributional assumptions made in DESeq and related tools. For example, in DESeq2, counts are modeled "as following a negative binomial distribution". This assumption is not valid when summing estimated counts of transcripts to obtain gene level counts, hence the justified concern of Michael Love that plugging in sums of estimated transcript counts could be problematic for DESeq2. In fact, even the estimated transcript counts themselves are not negative binomial distributed, and therefore also those are not appropriate for plugging into DESeq2. His concern is equally valid with many other "count based" differential expression tools.

>3. Fortunately there is a solution for performing valid statistical testing of differential abundance of individual transcripts, namely the method implemented in sleuth. The approach is described here. To test for differential abundance of genes, one must first address the question of what that means. E.g. is a gene differential if at least one isoform is? or if all the isoforms are? The tests of sleuth are performed at the granularity of transcripts, allowing for downstream analysis that can capture the varied questions that might make biological sense in specific contexts.

>In summary, please do not plug in rounded estimates of gene counts from kallisto into DESeq2 and other tools. While it is technically possible, it is not statistically advisable. Instead, you should use tools that make valid distributional assumptions about the estimates.

### Databases
* [ReCount is an online resource consisting of RNA-seq gene count datasets built using the raw data from 18 different studies](http://bowtie-bio.sourceforge.net/recount/)
* [The Digital Expression Explorer](http://dee.bakeridi.edu.au/index.html) The Digital Expression Explorer (DEE) is a repository of digital gene expression profiles mined from public RNA-seq data sets. These data are obtained from NCBI Short Read Archive.  
[blog post for it](http://genomespot.blogspot.com/2015/10/introducing-digital-expression-explorer.html?utm_content=bufferb6214&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)

### Gene Set enrichment analysis
* [Gene set analysis approaches for RNA-seq data: performance evaluation and application guideline](http://bib.oxfordjournals.org/content/early/2015/09/04/bib.bbv069.long)  
* 

### Pathway analysis
[Statistical analysis and visualization of functional profiles for gene and gene clusters: bioconductor clusterProfiler](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html) by GuangChuang Yu from University of HongKong. Can do many jobs and GSEA like figure. It is very useful and I will give it a try besides [GAGE](http://bioconductor.org/packages/release/bioc/html/gage.html).


### Fusion gene detection
* [fusioncatcher](https://github.com/ndaniel/fusioncatcher)  
* [PRADA](https://github.com/crazyhottommy/PRADA_pipeline_Verhaak_lab) from our lab

### Alternative splicing
* [SplicePlot: a tool for visualizing alternative splicing](http://montgomerylab.stanford.edu/spliceplot/index.html) Sashimi plots
* [Multivariate Analysis of Transcript Splicing (MATS)](http://rnaseq-mats.sourceforge.net/)
* [SNPlice](https://code.google.com/p/snplice/) is a software tool to find and evaluate the co-occurrence of single-nucleotide-polymorphisms (SNP) and altered splicing in next-gen mRNA sequence reads. SNPlice requires, as input: genome aligned reads, exon-intron-exon junctions, and SNPs. exon-intron-exon junctions and SNPs may be derived from the reads directly, using, for example, TopHat2 and samtools, or they may be derived from independent sources

### Allel specific expression
* paper [Tools and best practices for data processing in allelic expression analysis](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
* paper [Characterizing noise structure in single-cell RNA-seq distinguishes genuine from technical stochastic allelic expression](http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html)  
* paper [Tools and best practices for data processing in allelic expression analysis](http://www.genomebiology.com/2015/16/1/195)

### Single cell RNA-seq
* [On the widespread and critical impact of systematic bias and batch effects in single-cell RNA-Seq data](http://biorxiv.org/content/early/2015/08/25/025528) 
* [Ginkgo](http://qb.cshl.edu/ginkgo/?q=/ESjKTTeZIdnoGwEB4WTu) A web tool for analyzing single-cell sequencing data.
* [Seurat](http://www.satijalab.org/seurat.html) is an R package designed for the analysis and visualization of single cell RNA-seq data. It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
* [R package for the statistical assessment of cell state hierarchies from single-cell RNA-seq data](http://bioconductor.org/packages/devel/bioc/html/sincell.html)  
* [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) Differential expression and time-series analysis for single-cell RNA-Seq and qPCR experiments.
