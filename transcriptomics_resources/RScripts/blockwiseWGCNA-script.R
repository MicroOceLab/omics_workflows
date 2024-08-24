library(matrixStats)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(ggfortify)
library(WGCNA)
library(ggrepel)
library(patchwork)
library(dplyr)

# Check the each package for proper installation.
# Some can be installed via: 
#   install.packages("package name")
#       OR
#   BiocManager::install("package name")


##############
#
# Purpose of script: Run blockwise WGCNA pipeline
# Note: If you have a lot of features and computer memory becomes an issue,
#         this is an alternative approach to standard WGCNA. The author discusses
#         some details about the blockwise approach on the blog linked below:
#         https://peterlangfelder.com/2018/11/25/blockwise-network-analysis-of-large-data/
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
rv2 # Just to check the values of highest variances

# Plot ordered feature variance
# x = transcript number
# y = feature variance
theme_set(theme_bw(base_size = 10))
ggplot(rv2, aes(x=Seq,y=rowVars)) + geom_line() + scale_y_log10() +
  ggtitle("vst-transformed counts ordered by rowVar")

# Create new column containing variance
trans_filt_expr_df = as.data.frame(trans_filt_expr_mat)
trans_filt_expr_df$variance = rv

# Select only top 'n' genes with highest variance
# The value 0.25 in the code below says that only the top 25%
#   transcripts with highest variances will be selected. Change
#   this accordingly
varfilt_expr_df = trans_filt_expr_df[order(trans_filt_expr_df$variance, decreasing=TRUE), ][0:round(dim(trans_filt_expr_df)[1]*0.25, digits=0), ]
varfilt_expr_df = varfilt_expr_df[, !(names(varfilt_expr_df) %in% c("variance"))]

dim(varfilt_expr_df) # check dimension of filtered count matrix

# Transpose again
varfilt_expr_df = t(varfilt_expr_df)

# Uncomment the line below to save variance filtered data
#saveRDS(varfilt_expr_df, "1-trials/top25pct_genes_varfilt_expr_df.rds")

#Or simply load the variance filtered expression matrix
varfilt_expr_df = readRDS("1-trials/top25pct_genes_varfilt_expr_df.rds")


# ===== Choosing soft-thresholding (beta) parameter =====

# Check out WGCNA's paper and documentation to know more about the beta parameter:
#   Relocated manual: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0
#     The FAQ page gives some great advice! :)
#   Paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559


# Create a list of what beta values to test
powers = c(c(4:10), seq(from=12, to=30, by=2))

# Perform network topology analysis using different beta
sft = pickSoftThreshold(
  varfilt_expr_df, 
  powerVector=powers,
  networkType="signed",
  verbose=5)

sft

# Uncomment the line below if you want to save the variable 'sft'
#   so you won't have to re-run it next time - you can just import
#   the RData file instead.
# save(sft, file="1-trials/top_50pct_pickSoftThreshold_output.RData")

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


# ===== Block-wise network construction =====

# Construct network using blockwise approach
# Adjust maxBlockSize to a value that your memory can manage
# Adjust minModuleSize
# You can play around with different values of deepSplit (0 to 4) where
#   0 is the least sensitive and 4 is the most sensitive in
#   cutting tree (i.e. determining modules)
# For visualization purpose, it may be better to set numericLabels=FALSE, but
#   you may have to adjust some of the codes below since this was ran with TRUE

bwnet = blockwiseModules(
  varfilt_expr_df,
  maxBlockSize=20000,
  power=14,
  networkType="signed",
  TOMType="signed",
  saveTOMs=TRUE,
  saveTOMFileBase="bw-TOM",
  minModuleSize=20,
  reassignThreshold=0,
  deepSplit=0,
  pamRespectsDendro=FALSE,
  verbose=3,
  numericLabels=TRUE
)

# Uncomment the line below to save the result
#save(bwnet, file="1-trials/top25pct_genes_power16_blockwiseModules_output.RData")

#Check number of modules and members
table(bwnet$colors)


# ===== Plot dendogram =====

#Convert numerical labels to colors
merged_colors = labels2colors(bwnet$colors)
merged_colors

#plot dendogram
plotDendroAndColors(
  bwnet$dendrograms[[3]], 
  merged_colors[bwnet$blockGenes[[3]]],
  "Module colors",
  dendroLabels=FALSE,
  hang=0.03,
  addGuide=TRUE,
  guideHang=0.05)


# ===== Relate modules and MEs to external trait =====

# The codes below uses a different approach in finding modules 
#   associated with an external trait. Specifically, it uses linear models
#   to determine if the trait of interest is associated with a module.
#   It uses two models and the consensus is used to infer 
#   toxicity-associated modules.
# You can also refer to the WGCNA tutorial: Consensus-RelateModsToTraits.pdf
#   for a more direct approach

# Convert Toxic metadata column to numerical
# This may or may not be applicable in your metadata
Tox = ifelse(metadata$Toxic=="Y", 1, 0)
Tox

# Produce a model matrix that takes into account the following variables:
#   Toxic + Cross
#   Change these accordingly
mod_matr = as.data.frame(model.matrix(~ Toxic + Cross, metadata))

# MODEL 1:
# Create dataframes to store results
glm_mod_cor_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(glm_mod_cor_df) <- c("module", "pval", "coeff", "adj_rsq")

# Check association between module eigengenes (ME) and trait (toxicity) taking into account parental cross
# This loops through all module eigengenes
for (ME_num in names(bwnet$MEs)) {
  mod_trait_cor = lm(bwnet$MEs[, ME_num] ~ mod_matr$CrossT1xNT1 + mod_matr$CrossT1xT2 + mod_matr$CrossT2xNT2 + mod_matr$ToxicY)

  p_val = summary(mod_trait_cor)$coefficients["mod_matr$ToxicY", "Pr(>|t|)"]
  coeff = summary(mod_trait_cor)$coefficients["mod_matr$ToxicY", "Estimate"]
  adj_rsq = summary(mod_trait_cor)$adj.r.squared
  
  glm_mod_cor_df <- rbind(
    glm_mod_cor_df, 
    data.frame(module=ME_num, pval=p_val, coeff=coeff, adj_rsq=adj_rsq)
  )
  
  if (p_val <= 0.05) {
    print(paste(ME_num, p_val, coeff, adj_rsq))
  }
}

# Adjust p-value for multiple hypothesis testing
padj <- p.adjust(glm_mod_cor_df$pval, method="BH")
glm_mod_cor_df$padj <- padj

# Sort rows based on padj (adjusted p-value)
glm_mod_cor_df <- glm_mod_cor_df[order(glm_mod_cor_df$padj), ]
row.names(glm_mod_cor_df) <- NULL

# Get rows with significant p-adj value
glm_mod_cor_sig <- glm_mod_cor_df[glm_mod_cor_df$padj<=0.05, ]
glm_mod_cor_sig #Check modules significantly correlated with trait

# MODEL 2:
# Create dataframes to store results
lm_mod_cor_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(lm_mod_cor_df) <- c("module", "pval", "coeff", "adj.rsq")

# Check association between module eigengenes (ME) and trait (toxicity)
# This loops through all module eigengenes
for (ME_num in names(bwnet$MEs)) {
  mod_trait_cor <- lm(unlist(bwnet$MEs[, ME_num]) ~ mod_matr$ToxicY)
  
  p_val <- summary(mod_trait_cor)$coefficients["mod_matr$ToxicY", "Pr(>|t|)"]
  coeff <- summary(mod_trait_cor)$coefficients["mod_matr$ToxicY", "Estimate"]
  adj_rsq <- summary(mod_trait_cor)$adj.r.squared
  
  lm_mod_cor_df <- rbind(
    lm_mod_cor_df, 
    data.frame(module=ME_num, pval=p_val, coeff=coeff, adj.rsq=adj_rsq)
  )
  
  if (p_val <= 0.05) {
    print(paste(ME_num, p_val, coeff))
  }
}

# Adjust p-value for multiple hypothesis testing
padj <- p.adjust(lm_mod_cor_df$pval, method="BH")
lm_mod_cor_df$padj <- padj

# Sort rows based on padj (adjusted p-value)
lm_mod_cor_df <- lm_mod_cor_df[order(lm_mod_cor_df$padj), ]
row.names(lm_mod_cor_df) <- NULL

# Get rows with significant p-adj value
lm_mod_cor_sig <- lm_mod_cor_df[lm_mod_cor_df$padj<=0.05, ]
lm_mod_cor_sig #Check modules significantly correlated with trait

# Get modules significantly associated with trait in both models
trait_assoc_mods <- intersect(glm_mod_cor_sig$module, lm_mod_cor_sig$module)
trait_assoc_mods

# Hereafter, you can proceed with analyzing the functions in these modules,
#   or you can also do other statistical tests to determine certain functions
#   that are overrepresented in these modules relative to the transcriptome
#   e.g. overrepresentation analysis, fisher's exact test, etc.


# ===== Identify genes with high gene significance (GS), and module membership (MM) =====

# You may also want to explore certain genes that may have important 
#   biological roles in the modules. To do that you can look at gene 
#   significance (GS), and module membership (MM)
# GS is the correlation of feature and trait
# MM is the correlation of feature and module eigengene

# Correlate expression to MEs for MM
gene_mod_membership = as.data.frame(cor(varfilt_expr_df, bwnet$MEs, use="p"))

# Correlate expression to trait for GS
gene_trait_signif = as.data.frame(cor(varfilt_expr_df, Tox, use="p"))

#Just quick check for GOIs
gene_mod_membership["comp78043_c0_seq1", paste0("ME", MOI_num)]
gene_trait_signif["comp78043_c0_seq1", ]

# Combine MM and GS dataframes
# Select first which module you want to specifically analyze
MOI_num <- 66

gs_mm_df = data.frame(
  MM=gene_mod_membership[names(which(bwnet$colors==MOI_num)), paste0("ME", MOI_num), drop=FALSE],
  GS=gene_trait_signif[names(which(bwnet$colors==MOI_num)), "V1", drop=FALSE]
)
names(gs_mm_df) <- c("MM", "GS")

# MM vs GS plot
# Features on the upper right section (High MM and High GS) may
#   be interesting to further explore
gs_mm_df %>% ggplot(aes(x=MM, y=GS)) +
  geom_point(
    aes(alpha=0.1), 
    show.legend=FALSE,
    size=2
  ) +
  labs(title=paste0("MM vs GS: Module ", MOI_num)) +
  xlab("Module membership") +
  ylab("Gene significance") +
  theme_classic(base_size=16) +
  geom_text_repel(
    data=gs_mm_df, 
    aes(label=row.names(gs_mm_df)), 
    size=3,
    max.overlaps=17,
    position=position_dodge(width=0.5)
  )
  

gs_mm_df[DEGs, ]
verboseScatterplot(
  gene_mod_membership[row.names(GOI), paste0("ME", MOI_num)],
  gene_trait_signif[row.names(GOI), 1],
  xlab=paste0("module membership in ", MOI_color, " module"),
  ylab="gene significance for toxicity",
  main="module membership vs gene significance",
  cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=MOI_color)

text(
  gene_mod_membership[DEGs, paste0("ME", MOI_num)],
  gene_trait_signif[DEGs, 1],
  labels=DEGs,
  cex=0.7,
  pos=2)

abline(v=0.85, h=0.85, lty="dashed")

#Calculate intramodular connectivity
adj_mat = ((1-cor(varfilt_expr_df, use="p"))/2)^14
intramodularConnectivity()



