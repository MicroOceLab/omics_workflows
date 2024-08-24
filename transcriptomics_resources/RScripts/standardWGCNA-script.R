library(WGCNA)
library(matrixStats)
library(ggplot2)
library(pheatmap)

# Check the each package for proper installation.
# Some can be installed via: 
#   install.packages("package name")
#       OR
#   BiocManager::install("package name")


##############
#
# Purpose of script: Run standard WGCNA pipeline
# Note: WGCNA can be a memory intensive task if you have plenty of
#       features (transcripts, OTUs, ASVs, etc) in your data.
#       If you run into memory issues, you can check out the blockwise
#       WGCNA approach.
# Note: Replace the directories listed below as needed
#
##############

# ===== Set working directory =====

# All the relevant files should be placed relative to this directory
setwd("E:/Organized_Stuff/Work_Files/MSI/HABs_Watch/individual_research/R-processing/6-WGCNA")


# ===== Prep and import initial data =====

# Load expression matrix file
# This is the variance stabilized matrix generated from DESeq2. 
expr_mat = read.csv(
  "0-init_data/vst_counts.tsv", 
  sep="\t", 
  row.names=1
)
expr_mat = t(expr_mat)

dim(expr_mat) # samples x features

# Check if all samples and features are okay (i.e. zero variance,
#   missing entries, etc.)
# You will see from the printed output how many features WGCNA
#   deems to be 'bad'
gsg = goodSamplesGenes(expr_mat, verbose=3)
gsg$goodSamples #check if all samples are good

# Remove features that may not have information value for a
#    correlation analysis like WGCNA
filt_expr_mat = expr_mat[, gsg$goodGenes]
dim(filt_expr_mat)

# Cluster samples to check if there are outliers
# PCA described in DESeq2 script should also work for this
# If you find that there are outlier samples, you may omit
#   them to reduce noise in the data
sample_tree = hclust(dist(filt_expr_mat), method="average")
plot(sample_tree, main="Sample clustering to check for outliers")

# Load metadata
metadata = read.csv("0-init_data/samples_desc_full.txt", sep="\t", row.names=4)
metadata


# ===== (OPTIONAL) Filter expression matrix based on variance =====

# Removing low-variance features may help in (1) reducing
#   the dimension of the dataset, and thus, reducing memory
#   consumption as well, (2) removing transcripts that
#   do not have much information value (i.e. difficult to
#   correlate if counts across samples barely vary).
# There is no consensus on what threshold to use, so just
#   play aroudn and explore how much the resulting network
#   is affected by changing threshold.
# The first few codes below will help you visualize how
#   feature variance is distributed.

# Transpose count matrix
trans_filt_expr_mat = t(filt_expr_mat)
dim(trans_filt_expr_mat)

# Check distribution of variance
rv = rowVars(as.matrix(trans_filt_expr_mat))
rv2 = data.frame(Seq=seq(1:nrow(trans_filt_expr_mat)), rowVars=rv[order(rv, decreasing=TRUE)])
rv2

# Plot ordered feature variance
# x = transcript number
# y = feature variance
theme_set(theme_bw(base_size = 10))
ggplot(rv2, aes(x=Seq,y=rowVars)) + geom_line() + scale_y_log10() +
  ggtitle("vst-transformed counts ordered by variance")

# Create new column containing variance
trans_filt_expr_df = as.data.frame(trans_filt_expr_mat)
trans_filt_expr_df$variance = rv

# Select only top 'n' genes with highest variance
# Replace 50,000 by whatever value you want to try.
varfilt_expr_df = trans_filt_expr_df[order(trans_filt_expr_df$variance, decreasing=TRUE), ][0:50000, ]
varfilt_expr_df = varfilt_expr_df[, !(names(varfilt_expr_df) %in% c("variance"))]

# Transpose again
varfilt_expr_df = t(varfilt_expr_df)


# ===== Choosing soft-thresholding (beta) parameter =====

# Check out WGCNA's paper and documentation to know more about the beta parameter:
#   Relocated manual: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0
#     The FAQ page gives some great advice! :)
#   Paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559


# Create a list of what beta values to test
powers = c(c(4:10), seq(from=12, to=40, by=2))

# Perform network topology analysis using different beta
sft = pickSoftThreshold(
  varfilt_expr_df, 
  powerVector=powers, 
  verbose=5, 
  networkType="signed"
)

sft

# Uncomment the line below if you want to save the variable 'sft'
#   so you won't have to re-run it next time - you can just import
#   the RData file instead.
# save(sft, file="../R-processing/6-WGCNA/1-trials/varFilt_pickSoftThreshold_output.RData")

# Create plots to assess what beta parameter is suitable for your data
# Plot scale independence against beta
cex1 = 0.9

plot(
  sft$fitIndices[, 1], 
  -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
  xlab="Beta", ylab="Scale Free Topology Model Fit,signed R^2", 
  type="n",
  main=paste("Scale independence"))
text(
  sft$fitIndices[,1], 
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers, 
  cex=cex1, 
  col="red")

# Plot mean connectivity against beta
plot(
  sft$fitIndices[,1], 
  sft$fitIndices[,5], 
  xlab="Beta", ylab="Mean Connectivity", 
  type="n", 
  main=paste("Mean Connectivity"))
text(
  sft$fitIndices[,1], 
  sft$fitIndices[,5], 
  labels=powers, 
  cex=cex1,
  col="red")


# ===== Calc adjacency matrix and TOM =====

# Set soft_power based on selected beta
soft_power = 12

# Calculate adjacency matrix
# type='signed' takes into account the direction of correlation. Change
#   this value if this does not suit your intended objective
adj_mat = adjacency(
  varfilt_expr_df, 
  power=soft_power,
  type="signed")

# Calculate TOM
TOM = TOMsimilarity(adj_mat, TOMType="signed")

# Convert to dissimilarity matrix
dissTom = 1 - TOM
dim(dissTom) # Just checking dimensions


# ===== Cluster features into modules =====

# Create tree
gene_tree = hclust(as.dist(dissTom), method="average")

# View tree
plot(
  gene_tree, 
  xlab="", 
  sub="", 
  main="Gene clustering based on TOM-dissimilarity", 
  labels=FALSE)

# Set the minimum size of a module
min_mod_size = 20

# Cluster features into modules using a dynamic approach
dynamic_mods_ds0_pamstageF = cutreeDynamic(
  dendro=gene_tree,
  distM=dissTom,
  deepSplit=0, 
  pamRespectsDendro=FALSE,
  minClusterSize=min_mod_size,
  pamStage=FALSE
)

table(dynamic_mods_ds0_pamstageF)

dynamic_colors_ds0_pamstageF = labels2colors(dynamic_mods_ds0_pamstageF)
table(dynamic_colors_ds0_pamstageF)

plotDendroAndColors(
  gene_tree, dynamic_colors_ds0, "Dynamic Tree Cut",
  dendroLabels=FALSE,
  hang=0,
  addGuide=FALSE,
  guideHang=0.05,
  main="Gene dendogram and module colors"
)

# After generating the modules, you can refer to blockwiseWGCNA-script.R
#   for subsequent steps you can do to explore your data (e.g. module-trait correlation,
#   module membership, and gene significance)