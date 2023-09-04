# CNVmap
 CNVmap is a tool for visualizing Copy Number Variations (CNVs) at genome, chromosome and gene levels. The inputs of this tool are .cns and .cnr file formats obtained from the tool CNVkit (https://cnvkit.readthedocs.io/en/stable/).
 
CNVs, also called Copy Number Alterations (CNAs) are a subtype of Structural Variants (SVs). Althought their definition is not clearly determined, they are portions of the genome that can either be lost (deletions) or duplicated (amplifications).<br />
CNVs are a source of genetic variability, however they are also involve in development of some diseases. In particular, Somatic Copy Number Alterations contribute to cancer initiation, progression and therapeutic resistance. For example, the level of copy number alteration (CNA) is associated with recurrence of primary prostate cancer or with disease-free and overall survival in primary breast, endometrial, renal clear cell, thyroid, and colorectal cancer [1]. For these reasons, the detection and study of Copy Number Alterations could play a key role in the treatment of several tumors. **CNVkit** is a software which detect Copy Number Variants and Alterations genome wide [2]. In order to inspect CNAs at different genomic levels and produce publication-graphic figures we developed the R-package **CNVmap**. 

&nbsp;

## Installation

You will need to install it via _GitHub_ using the `devtools` package.<br /> 
To install `devtools`:

```r
install.packages("devtools")
library(devtools)
```

To install CNVmap and its dependencies:
```r
devtools::install_github("BrandoPoggiali/CNVmap")
```
&nbsp;

## Update
To update the package, use the following command:
```r
devtools::update("CNVmap")
```
&nbsp;

## Usage
Upload .cnr and cns. file (generated by CNVkit software) in R using the `data.table` library.

```r
library(data.table)

cnr_data <- fread("path/name.cnr", sep="\t")
cns_data <- fread("path/name.cns", sep="\t")
```
&nbsp;
### All chromosomes plot
The function for plotting CNAs in all chromosomes is time-intensive, it could take some minutes:

```r
plot_all_chr(cnr_data, cns_data, only_autosomal = TRUE, chr_text_size=5)
```

![](img/All_autosomes.png)


&nbsp;

### Single chromosome plot
To plot CNAs at chromosome level: 

```r
plot_single_chr(cnr_data, cns_data, chr = "chr6")
```
![](img/Chromosome_6.png)

&nbsp;

If you want to display locations of a set of genes of interest in the selected chromsome: 

```r
genes_list <- c("BRCA2","FLT1","NOL7", "TP53", "MYC", "NUPL1", "POMP", "HLA-A", "SOX21", "ARG1", "MYO6", "ELOVL5")
plot_single_chr(cnr_data, cns_data, chr = "chr6", genes = genes_list, gene_text_size = 3.3)
```

![](img/Chromosome_6_genes.png)

&nbsp;

To add a chromosome icon showing deleted and amplified region:

```r
plot_single_chr(cnr_data, cns_data, chr = "chr6", genes = genes_list, gene_text_size = 3.3, chr_picture = TRUE)
```

![](img/Chromosome_6_genes_icon.png)

&nbsp;

### Single gene plot
It is possible to explore CNVs at gene level. You have to specify the name of a gene with the HGNC ID (the function automatically plot the canonical transcript according to Ensembl 107).

```r
plot_single_gene(cnr_data, cns_data, gene = "TP53")
```
![](img/TP53.png)

&nbsp;

The function allow to plot a detailed annotation with both all transcripts of the gene of interest and the regulatory elements present in that region. 
The color of the regulatory elements follow the in UCSC browser (red for promoter, orange for distal enhancer, yellow for proximal enhancer, pink for DNase-H3K4me3, and blue for CTCF-only).

```r
plot_single_gene(cnr_data, cns_data, gene = "TP53", all.transcripts = TRUE, regulatory.elements = TRUE)
```

![](img/TP53_annotated.png)

&nbsp;

## Citation 
If you use this tool, please consider citing our publication.
B. Poggiali *et al*(2023). Multiomic analysis of HER2-enriched and AR-positive breast carcinoma with apocrine differentiation and an oligometastatic course: a case report. *Frontiers in Oncology*.

&nbsp;

## Contact

Questions, suggestions, and bug reports are welcome and appreciated.
- **Author**: brando Poggiali
- **email**: brando.poggiali@iit.it
- **2° mail**: poggialibrando1995@gmail.com

&nbsp;

## References
[1] Hieronymus, Haley, et al. "Tumor copy number alteration burden is a pan-cancer prognostic factor associated with recurrence and death." Elife 7 (2018): e37294.<br />
<br />
[2]Talevich, Eric, et al. "CNVkit: genome-wide copy number detection and visualization from targeted DNA sequencing." PLoS computational biology 12.4 (2016): e1004873.
