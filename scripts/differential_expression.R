############################################################
# Differential Gene Expression Analysis using DESeq2
# Dataset: Wang Expression Set
# Comparison: Brain vs Other Tissue Types
############################################################

# Load dataset


download.file(
  "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/wang_eset.RData",
  destfile = "wang_eset.RData"
)

load("wang_eset.RData")

# ----------------------------------------------------------
# Load required packages
# ----------------------------------------------------------

library(Biobase)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ----------------------------------------------------------
# Prepare count matrix and sample metadata
# ----------------------------------------------------------

count_matrix <- exprs(wang.eset)[, 10:21]
col_data <- pData(wang.eset)[10:21, ]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ cell.type
)

# ----------------------------------------------------------
# Create Brain vs Other comparison
# ----------------------------------------------------------

dds$type <- ifelse(
  dds$cell.type %in% c("cerebellum", "mixed.brain"),
  "brain",
  "other"
)

dds$type <- factor(dds$type)
dds$type <- relevel(dds$type, ref = "other")

design(dds) <- ~ type

# Remove low-count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# ----------------------------------------------------------
# Run DESeq2 analysis
# ----------------------------------------------------------

dds <- DESeq(dds)

results_default <- results(dds)

summary(results_default)

# ----------------------------------------------------------
# Visualize results
# ----------------------------------------------------------

plotMA(results_default)

plotMA(
  results_default,
  ylim = c(-10, 10)
)

# ----------------------------------------------------------
# Most significant gene
# ----------------------------------------------------------

top_gene <- rownames(results_default)[which.min(results_default$padj)]

top_gene

# ----------------------------------------------------------
# Differential expression with effect size threshold
# ----------------------------------------------------------

results_lfc <- results(
  dds,
  lfcThreshold = 2,
  alpha = 0.10
)

results_lfc

table(results_lfc$padj < 0.10)

plotMA(
  results_lfc,
  ylim = c(-10, 10)
)

summary(results_lfc)

sum(
  results_lfc$padj < 0.10 &
    results_lfc$log2FoldChange > 2,
  na.rm = TRUE
)

# ----------------------------------------------------------
# Top 20 genes ranked by Wald statistic
# ----------------------------------------------------------

top_genes <- rownames(results_default)[
  head(order(results_default$stat, decreasing = TRUE), 20)
]

top_gene_ensembl <- top_genes[1]
top_gene_ensembl

# ----------------------------------------------------------
# Annotate top gene
# ----------------------------------------------------------

top_gene_name <- mapIds(
  org.Hs.eg.db,
  keys = top_gene_ensembl,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top_gene_name
