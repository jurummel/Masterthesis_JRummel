# Masterthesis Julian Rummel
## Title
Combining epigenomic and transcriptomic data for the analysis of allele-specific regulation in a heart failure model

## Description
This is a GitHub repository to collect all results with the corresponding code to make it accessible and reproducible for readers.

## File Hierarchy 
- Results
  - mRNA-seq 
    contains all mRNA-seq results (Quality Control files, Differentially Expressed Genes (DEG), Clustering, Functional Enrichment Analysis)
  - ATAC-seq
    contains all ATAC-seq results (Quality Control Files, Peaks, Differential Accessibility Region (DAR))
  - Combined
    contains all results from STARE and HOMER (Gene-Enhancer Interactions, Motif Enrichment Analysis)
  
- Code
  - contains code for Gene Expression Analysis with DESeq2, Clustering of DEGs, Functional Enrichment Analysis, and Analysis of DARs
   
- Thesis
  - Contains the PDF-file of the thesis

## Dependencies and Sources
The analyses depend on following R libraries (all versions are found in the sessionInfo): 
- [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html)
- [VennDetail](https://www.bioconductor.org/packages/release/bioc/html/VennDetail.html)
- [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
- [apeglm](https://bioconductor.org/packages/release/bioc/html/apeglm.html)
- [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html)
- [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
- [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
- [dendextend](https://cran.r-project.org/web/packages/dendextend/index.html)
- [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2/index.html)
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)

The analyses is also dependend on the following tools:
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) (version 22.9.0)
- [snakePipes](https://snakepipes.readthedocs.io/en/latest/content/setting_up.html) (version 2.1.1)
- [samtools](http://www.htslib.org) (version 1.6)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 1.2.3)
- [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (version 0.6.7)
- [macs2](https://docs.csc.fi/apps/macs2/) (version 2.2.6)
- [homer](http://homer.ucsd.edu/homer/motif/) (version 4.11)
- [stare](https://stare.readthedocs.io/en/latest/Main.html) (version 1.0.3)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.9)
- [MultiQC](https://multiqc.info) (version 1.0.dev0)

The code is inspired by the following tutorials:
- [Analyzing RNA-seq data with DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) (Bioconductor, authors: Michael I. Love, Simon Anders, and Wolfgang Huber)
- [Exploring gene expression patterns using clustering methods](https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html) (Github Pages, authors: Hugo Tavares and Georg Zeller)
- [Differential Peak calling using DiffBind](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html) (GitHub Pages, authors: Meeta Mistry - hbctraining)

## Author
  - Julian Rummel
