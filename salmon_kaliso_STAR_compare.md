### Why I am testing again?

I know there are papers comparing different RNA-seq pipelines. For example:
[A benchmark for RNA-seq quantification pipelines](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0940-1)

I just got several RNA-seq data to play with, and I think it is a good time-point for me to get my hands wet on those RNA-seq quantification tools (especially those alignment-free ones) and get a personal idea of how different tools perform. I am not 
doing bench-marking, as one should simulate the RNA-seq reads by e.g. [polyester](https://github.com/alyssafrazee/polyester) to have 
the ground truth.

Choosing alignment based tools (such as tophat, STAR, bowtie, HISAT) or alignment free ones depends on the purpose of your study. `Salmon` and `kallisto` requires the reads "pesudo-map" to the transcriptome, so one has to provide a fasta file containing all the transcripts you want to quantify. Therefore, if you want to find novel transcripts, you probably should go with the alignment based methods. It is also shown recently that [Widespread intron retention diversifies most cancer transcriptomes](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0168-9). If you want to do similar things, you need to use mappers with the genome (not the transcriptiome) as a reference. 

`kallisto` can output a pseudo-bam which can be useful for some people. `Salmon` will have the same functionality in the next release according to Rob.

Now, let's begin my analysis.

### Testing [Salmon](https://github.com/COMBINE-lab/salmon) for RNA-seq quantification

Download the binary(v0.6.0) by:

```bash
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz
tar xvzf SalmonBeta-0.6.0_DebianSqueeze.tar.gz
```

Read the documentation on how to use it [here](http://salmon.readthedocs.io/en/latest/salmon.html#using-salmon)

To use Salmon in quasi-mapping-based mode, then you first have to build an Salmon index for your transcriptome.
[Ensemble](http://useast.ensembl.org/info/data/ftp/index.html) release 75 is the latest version for GRCh37 (hg19).
starting from release 76 it is GRCh38 (hg38)

Download the hg19 version of cDNA and non-coding RNA fasta:

```bash
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz  
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.75.ncrna.fa.gz
## merge together 

gunzip -c Homo_sapiens.GRCh37.75.cdna.all.fa.gz Homo_sapiens.GRCh37.75.ncrna.fa.gz > Homo_sapiens.GRCh37.75.cdna.ncrna.fa
```
First, you run the Salmon indexer:

>Default index --- The quasi index has been made the default type. This means that it is no longer necessary to provide the 
--type option to the index command. The fmd index remains enabled, but may be removed in a future version.

```bash
salmon index -t Homo_sapiens.GRCh37.75.cdna.ncrna.fa -i Homo_sapiens.GRCh37.75_quasi_index 
```
It finished in a little over 10 mins.

quantify the transcript:

library is single-end stranded, reads from reverse strand

```bash
salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat 3R_S18_L002_R1_001.fastq.gz) -o 3R_transcripts_quant

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat 50R_S19_L002_R1_001.fastq.gz) -o 50R_transcripts_quant

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat WT_S17_L002_R1_001.fastq) -o WT_transcripts_quant
```
The quantification finishes within minutes!

### quantification in gene level
Salmon can also give gene-level quantification as long as feed a gtf file 

```bash
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat 3R_S18_L002_R1_001.fastq.gz) -o 3R_transcripts_quant -g ~/annotations/Homo_sapiens.GRCh37.75.gtf

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat 50R_S19_L002_R1_001.fastq.gz) -o 50R_transcripts_quant -g ~/annotations/Homo_sapiens.GRCh37.75.gtf

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat WT_S17_L002_R1_001.fastq.gz) -o WT_transcripts_quant -g ~/annotations/Homo_sapiens.GRCh37.75.gtf
```
It is recommended using [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html) to get the gene-level quantification. I asked the difference between `tximport` and `salmon quant -g`.
The developer of salmon @[Rob Patro](https://twitter.com/nomad421) answered:
>Main diffs I can think of (1) in R   
(2) integrated with DESeq2  
(3) Can derive multi-sample effective gene lengths

### testing [kallisto](https://pachterlab.github.io/kallisto/starting) for quantification.
Build index first:

```bash
kallisto index -i Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx Homo_sapiens.GRCh37.75.cdna.ncrna.fa
```
It took me around 15 mins to build the index.

quantification:

`-l, --fragment-length=DOUBLE  Estimated average fragment length`  
`-s, --sd=DOUBLE               Estimated standard deviation of fragment length`
                              `(default: value is estimated from the input data)`
>In the case of single-end reads, the -l option must be used to specify the average fragment length. Typical Illumina libraries produce fragment lengths ranging from 180–200 bp but it’s best to determine this from a library quantification with an instrument such as an Agilent Bioanalyzer.

see a question by [James](https://twitter.com/JamesPBLloyd) [on the google group](https://groups.google.com/forum/#!searchin/kallisto-sleuth-users/single$20end/kallisto-sleuth-users/VPJfzL502bw/e2JDq7ezBgAJ)

>Common values for single end reads are insert length 200 and sd 20. If you have any better information, like the person who prepped the library or better yet, data from bioanalyzer that will of course be better. 

`kallisto` can take `.gz` files.

```bash

kallisto quant -t 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o 3R_kaliso_output --single -l 200 -s 20 3R_S18_L002_R1_001.fastq.gz

kallisto quant -t 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o 50R_kaliso_output --single -l 200 -s 20 50R_S19_L002_R1_001.fastq.gz

kallisto quant -t 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o WT_kaliso_output --single -l 200 -s 20 WT_S17_L002_R1_001.fastq.gz
```
Finished in ~6 mins. again, blazing fast as `Salmon` does.

### compare with Kalisto and STAR-HTseq
