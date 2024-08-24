library(DESeq2)
library(tximport)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(gridExtra)
library(reshape2)
library(hash)
library(readr)
library(ggplot2)

# Check the each package for proper installation.
# Some can be installed via: 
#   install.packages("package name")
#       OR
#   BiocManager::install("package name")


##############
#
# Purpose of script: Run DESeq2 to perform differential expression analysis
# Note: This script was originally used to infer genes differentially expressed
#       between toxic and non-toxic strains of A. minutum. Change the parameters
#       of the functions below accordingly to make it fit your objectives
# Note: Replace the directories listed below as needed
#
##############


# ===== Set working directory =====

# All the relevant files should be placed relative to this directory
setwd("E:/Organized_Stuff/Work_Files/MSI/HABs_Watch/individual_research")


# ===== Prep and import initial data =====

# If you want to analyze at gene level, uncomment the code below, and
#   supply 'gene2trans' to the 'tx2gene' parameter of the 'tximport' function below.
#   Otherwise, proceed with the script as it is
# gene2trans = read.csv('./R-processing/gene2trans.tabular', header=TRUE, sep='\t') 

# Load metadata; TSV file with 1st row as column names
metadata_full <- read.csv(
  './R-processing/samples_desc_full.txt', 
  header=TRUE, 
  sep='\t'
)

# Set filepath for RSEM outputs
# You may need to modify the code below to set the 
#   proper file path of RSEM outputs
files_full <- file.path(
  './data/temp_wsl_storage/initial_exploration/4-mapping',
  paste(
    substring(metadata_full$Strain, 1, 30),
    '/rsem/RSEM.isoforms.results', 
    sep=''
  )
)

# Assign sample names to file paths
names(files_full) <- metadata_full$Strain

# Import RSEM output files
txi_full <- tximport(
  files_full, 
  type='rsem', 
  txOut=TRUE
)

# Import data as DESeq2 dataset
# The 'design' argument dictates what covariates (metadata variables)
#   you have that may be affecting the expression of transcripts.
#   You can do prior analysis like PCA (see below) to decide which variables
#   to include, or you can base it on previous knowledge/studies.
dds_full <- DESeqDataSetFromTximport(
  txi=txi_full,
  colData=metadata_full,
  design=~Toxic+Cross
)

# Check dimensions. It should be t x s
#   where t is the number of transcripts, and
#   s is the number of samples
dim(counts(dds_full))


# ===== Run DESeq2 analysis pipeline =====

# Run DESeq2
dds_full <- DESeq(dds_full)

# Plot dispersion estimates
# It's often used for diagnostic purposes, know more in the link below
#   https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#dispersion-plot-and-fitting-alternatives
plotDispEsts(dds_full)

# Check what possible pairwise comparisons
resultsNames(dds_full)

# Extract DESeq2 results
# Adjust lfcThreshold (log fold change) and alpha (p-value) accordingly
res_full <- results(
  dds_full, 
  lfcThreshold=1, 
  alpha=0.05, 
  name="Toxic_Y_vs_N"
)

# Print full results table
res_full

# Summarize results
summary(res_full)

# Get results table for DEGs only
nonNA_res_full <- res_full[!is.na(res_full$log2FoldChange) & !is.na(res_full$padj), ]
DEGs_res_full <- nonNA_res_full[abs(nonNA_res_full$log2FoldChange) > 1 & nonNA_res_full$padj < 0.05, ]

DEGs_res_full

# Get names of DE transcripts
DEG_names_full <- row.names(DEGs_res_full)
DEG_names_full

# Check model-fitted coefficients
DEGs_coef <- coef(dds_full)[DEG_names_full,]
DEGs_coef

# After running this section, you can check out the 2 bottomost sections for 
#   some ways to explore your data


# ===== OPTIONAL (READ DESCRIPTION FIRST): Shrink LFCs =====

# If you find that a lot of transcripts are differentially expressed,
#   it might be helpful to use lfcShrink to shrink fold changes of transcripts that may
#   be inaccurate or has higher noise (partly due to low expression values). Otherwise, you may
#   skip this entire section
# If you decide to run this section, make sure to change the relevant variables
#   in the subsequent sections by the new variables generated in this section

# Shrink LFC
res_shrunk_full <- lfcShrink(dds_full, coef="Toxic_Y_vs_N", type="apeglm", res=res_full)

# Get DEGs LFC-shrunk results
nonNA_res_shrunk_full <- res_shrunk_full[!is.na(res_shrunk_full$log2FoldChange), ]
DEGs_res_shrunk_full <- nonNA_res_shrunk_full[abs(nonNA_res_shrunk_full$log2FoldChange) > 1 & nonNA_res_shrunk_full$padj < 0.05, ]

DEGs_res_shrunk_full

# Extract DE transcripts
DEG_names_shrunk_full <- row.names(DEGs_res_shrunk_full)
DEG_names_shrunk_full

#Check coefficients of DEGs
DEGs_shrunk_coef <- coef(dds_full)[DEG_names_shrunk_full,]


# ==== Plot DEG heatmap =====

# Extract normalized count matrix
vsd_mat_full <- assay(vsd_full)

# Calculate correlation matrix
vsd_cor_full <- cor(vsd_mat_full)

# Plot DEG heatmap with metadata labels
pheatmap(
  vsd_mat_full[row.names(DEGs_res_full), ], 
  main="Heatmap of expression levels (vst) of DEGs", 
  cluster_rows=FALSE, 
  cluster_cols=TRUE,
  show_rownames=FALSE,
  annotation_col=as.data.frame(colData(dds_full)[, c("Cross", "Toxic")])
)


# ===== Boxplot for DEGs =====

# Subset the table to only include the DEGs
# You may need to subset further if you have a lot of DEGs
#   because the resulting plot may be too crowded
gene_expr_full <- vsd_mat_full[unlist(row.names(DEGs_res_full)), ]

# Melt table and rename columns
gene_expr_melt_full <- melt(gene_expr_full, na.rm=TRUE)
colnames(gene_expr_melt_full)[1] <- "DEGs"
colnames(gene_expr_melt_full)[2] <- "Strain"

gene_expr_melt_full

# Attach metadata information
gene_data_full <- merge(gene_expr_melt_full, metadata_full, by.x="Strain", by.y="Strain")
gene_data_full

# Plot boxplots
ggplot_boxplot_full <- ggplot(gene_data_full, aes(x=Toxic, y=value)) +
  geom_boxplot(outlier.shape=NA) +
  labs(title=paste("Boxplots of DEG expression levels")) +
  geom_point(aes(color=Cross), position=position_dodge(width=0.5), size=2) +
  facet_wrap(~DEGs, scales="free")

ggplot_boxplot_full


# ===== Generate MA Plot =====

# Produce MA plot (normalized counts vs LFC). DEGs found at the upper and 
#   lower right (i.e. high counts and high abs(LFC)) typically convey
#   stronger signals
# It somtimes helps if you use the LFC-shrunk table here so that
#   the you see less noise in the graph
# See the link below for reference:
#   https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results
plotMA(res_shrunk_full, ylim=c(-2, 2), main="MA plot", alpha=0.05)

# Custom MA plot using ggplot2
# In here, let's include in the visualization the LFC standard error (lfcSE).
#   Ideally, the DEGs should have smaller lfcSE (smaller triangles) which
#   indicates more narrow confidence interval
res_full_df <- as.data.frame(res_full)
res_full_df <- na.omit(res_full_df)

res_full_df$lfc

ggplot(res_full_df, aes(baseMean, log2FoldChange, size=lfcSE)) +
  geom_point(
    na.rm=TRUE,
    color=ifelse(res_full_df$padj < 0.05, "red", "grey20"),
    alpha=0.5,
    shape=ifelse(abs(res_full_df$log2FoldChange) > 1, 17, 16)) +
  scale_x_log10() +
  xlab("Mean of normalized counts") +
  ylab("log fold change") +
  labs(title="MA plot") +
  theme_minimal() +
  scale_size_area(limits=c(min(res_full_df$lfcSE), max(res_full_df$lfcSE)))


# ===== Generate Volcano Plot =====

# Produce volcano plots: LFC vs log(p-value). This is neat and common way of summarizing
#   the results of DE analysis

# For unlabeled volcano plot
# Adjust the parameters accordingly, some important arguments to change:
#   pCutoff = p-value cutoff
#   FCcutoff = fold-change cutoff
# Check the other arguments in their documentation or through help option
EnhancedVolcano(
  res_full,
  lab=NA,
  x='log2FoldChange',
  y='padj',
  pCutoff=5e-02,
  FCcutoff=1,
  labSize=2.5,
  title="Volcano Plot",
  subtitle=NULL,
  drawConnectors=TRUE,
  arrowheads=FALSE,
)

# For labeled volcano plot, do the following steps:

# 1. Make a copy of the results table so that we don't modify the original one
res_full_copy <- res_full

# 2. For downregulated transcripts, rename transcript IDs by gene name
# In the gsub function, the first argument is the transcript ID while the
#   second argument is the gene name. It's up to you which data points you
#   want labeled
rownames(res_full_copy) <- gsub("comp120540_c0_seq1", "CspA", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp63732_c0_seq2", "Ubiquitin C", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp116212_c0_seq3", "EF-hand", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp109639_c0_seq2", "Cysteine desulfurase", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp89316_c0_seq2", "LSU ribosomal protein L17e", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp97267_c1_seq1", "MFS family permease", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp115950_c0_seq2", "RND superfamily", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp110574_c0_seq2", "CobT", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp68927_c4_seq7", "Ubiquitin C", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp119614_c0_seq2", "Chk2", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp119614_c0_seq4", "Chk2", row.names(res_full_copy))

# 2. For top 10 upregulated transcripts, rename transcript IDs by gene name
# In the gsub function, the first argument is the transcript ID while the
#   second argument is the gene name. It's up to you which data points you
#   want labeled
rownames(res_full_copy) <- gsub("comp106437_c1_seq3", "Inner membrane complex protein", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp118540_c0_seq28", "Chromosome segregation ATPase", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp119515_c0_seq6", "NIMA-related kinase", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp125505_c0_seq4", "BamB", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp98170_c0_seq1", "Tubulin folding cofactor D", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp75290_c0_seq3", "Dip2/Utp12 Family", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp98774_c0_seq1", "WD40 repeat", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp84842_c0_seq2", "SET domain-containing protein", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp102365_c0_seq2", "Degradation arginine-rich protein for mis-folding", row.names(res_full_copy))
rownames(res_full_copy) <- gsub("comp79259_c0_seq1", "Poly-adenylate binding protein, unique domain", row.names(res_full_copy))

# 3. Create a list of labels you want to include. In here, those that contain
#   the string 'comp' will not be included. This leaves us with the rownames
#   with gene labels
volc_plot_labels <- unique(rownames(res_full_copy)[!grepl("^comp", rownames(res_full_copy))])

# 4. Generate the labeled volcano plot
EnhancedVolcano(
  res_full_copy,
  lab=rownames(res_full_copy),
  selectLab=volc_plot_labels,
  max.overlaps=Inf,
  boxedLabels=TRUE,
  labSize=4,
  x='log2FoldChange',
  y='padj',
  pCutoff=5e-02,
  FCcutoff=1,
  title="Volcano Plot",
  subtitle=NULL,
  drawConnectors=TRUE,
  widthConnectors=1.5,
  colConnectors='black',
  arrowheads=FALSE,
  raster=TRUE
)


# ===== Plot PCA (for data exploration) =====

# Normalize table using variance stabilizing transformation
vsd_full <- vst(dds_full, blind=TRUE)

# Make instances of plotPCA which can be used for ggplot modifications
pca_full <- plotPCA(
  vsd_full, 
  intgroup=c('TPG', 'Cross', 'Toxic', 'Strain'), 
  returnData=TRUE
)

#Create ggplot PCA plots
pca_plot_full <- ggplot(pca_full, aes(PC1, PC2, color=Toxic, shape=Cross, label=Strain)) +
  ggtitle("Full") +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100 * attr(pca_full, 'percentVar'))[1],"% variance")) +
  ylab(paste0("PC2: ",round(100 * attr(pca_full, 'percentVar'))[2],"% variance")) +
  coord_fixed() +
  geom_text_repel(size=3, colour="black", max.overlaps=50)

pca_plot_full


# ===== Plot sample correlation heatmap and DEG expression heatmap (for data exploration) =====

# Normalize table using variance stabilizing transformation
vsd_full <- vst(dds_full, blind=TRUE)

# Extract normalized count matrix
vsd_mat_full <- assay(vsd_full)

# Export the variance stabilized count matrix
# This can be used as input to WGCNA
write.table(vsd_mat_full, "R-processing/5-DE_analysis/Full_shrunken_data/vst_counts.tsv", sep="\t", quote=FALSE)

# Calculate correlation matrix
vsd_cor_full <- cor(vsd_mat_full)

# Plot pairwise sample correlation heatmap
pheatmap(vsd_cor_full, main='Sample correlation heatmap')
