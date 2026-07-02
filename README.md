RNA-Seq Differential Expression Analysis
This process performs differential gene expression analysis using DESeq2 on the Wang Expression Set dataset

PROJECT OVERVIEW
The goal is to identify differentially expressed genes between brain tissues and other tissue types. 
Methods used: Data processing, Count filtering, DESeq2 normalization, Differential expression testing, Log2 fold-change thresholding, Gene annotation using org.Hs.eg.db, and MA plot visualization

How to Run
Install required packages: install.packages("BiocManager")
BiocManager::install(c("DESeq2", "Biobase","org.Hs.eg.db", "AnnotationDbi"))
Run the script: source("scripts/differential_expression.R")

KEY OUTPUT
Differential expression results, significant genes(adjusted p-value< 0.1), MA plot, Annotation of the top-ranked gene

Skills Demonstrated
RNA-seq analysis, statistical genomics. DESeq2 workflow, data visualization, Bioconductor annotation, and reproducible research practices
