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

salmon quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75_quasi_index -l SR -r <(zcat WT_S17_L002_R1_001.fastq) -o WT_transcripts_quant -g ~/annotations/Homo_sapiens.GRCh37.75.gtf
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

```bash

kallisto quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o 3R_kaliso_output --single -l 200 -s 20 3R_S18_L002_R1_001.fastq.gz

kallisto quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o 50R_kaliso_output --single -l 200 -s 20 50R_S19_L002_R1_001.fastq.gz

kallisto quant -p 10 -i ~/annotations/Homo_sapiens.GRCh37.75.cdna.ncrna.kalisto.idx -o WT_kaliso_output --single -l 200 -s 20 WT_S17_L002_R1_001.fastq
```
### compare with Kalisto and STAR-HTseq
