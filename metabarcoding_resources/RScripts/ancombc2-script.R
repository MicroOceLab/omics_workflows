library(ANCOMBC)
library(qiime2R)
library(phyloseq)
library(dplyr)
library(tidyr)
library(microbiome)
library(plyr)

# Check the each package for proper installation.
# Some can be installed via: 
#   install.packages("package name")
#       OR
#   BiocManager::install("package name")


##############
#
# Purpose of script: Run ANCOM-BC2 to infer differentially abundant features (ASVs/OTUs) or taxa
# Note: ANCOM-BC2 may take a while to execute
# Required Files:
#       (1) Feature table from QIIME2 in .qza format
#       (2) Tab-separated metadata file
#       (3) Tab-separated taxonomy table exported from QIIME2 FeatureData[Taxonomy]
# Reference tutorial:
#       https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#
##############

# ===== Set working directory =====

setwd("E:/Organized_Stuff/Work_Files/MSI/RB and DJP/Singapore data_Brian")


# ===== Load and process feature table =====

# Load qza file
# Replace filename/path appropriately
feature_data <- read_qza('qza/filter_table_result.qza')
feature_table <- feature_data$data


# ===== Load and process metadata =====

# Load metadata 
metadata <- read.table(
  'metadata/metadata_complete_UTF8.tsv', 
  sep='\t',
  header=1,
)

# Remove row 1 which just displays column data type (e.g. numerical, categorical, etc)
# Omit line of code below if irrelevant in your metadata
metadata <- metadata[-1, ] 

# Set sample.id as dataframe row.names
row.names(metadata) <- metadata$sample.id
metadata <- metadata[, -1] # remove sample.id which are now row names


# ===== Load and process taxonomy table =====

# Load feature taxonomy table
taxonomy_table <- read.table(
  'taxonomy/taxonomy.tsv',
  sep='\t',
  header=1
)

# Remove column 3: confidence column because it is irrelevant here
taxonomy_table <- taxonomy_table[, -3]

# Separate taxonomy levels into different columns
taxonomy_table <- taxonomy_table %>%
  separate(Taxon, into = c("Lvl1", "Lvl2", "Lvl3", "Lvl4", "Lvl5", "Lvl6", "Lvl7"), sep = "; ", fill = "right")

# Set feature IDs as row.names
row.names(taxonomy_table) <- taxonomy_table$Feature.ID
taxonomy_table <- taxonomy_table[, -1]
taxonomy_table <- taxonomy_table[row.names(feature_table), ]


# ===== Convert to tables to matrix =====

feature_table <- as.matrix(feature_table)
taxonomy_table <- as.matrix(taxonomy_table)


# ===== Convert to phyloseq object =====

feature_phyloseq <- otu_table(feature_table, taxa_are_rows=TRUE)
taxonomy_phyloseq <- tax_table(taxonomy_table)
metadata_phyloseq <- sample_data(metadata)

data_phyloseq <- phyloseq(feature_phyloseq, taxonomy_phyloseq, metadata_phyloseq)


# ===== Run ANCOMBC2 =====

# Set seed so that experiment is repeatable
set.seed(123)

# Run ANCOM-BC2 analysis
# The argument 'fix_formula' tells ANCOM-BC2 which metadata variables to 
#   include in the model. In this case, we're saying that the model should
#   take into account the type of substrate and timepoint in fitting ASV
#   counts. Change this to whatever applies in your case
# The 'group' argument specifies which variable to perform the comparisons.
#   This must be a categorical variable. Again, change this as needed.
output = ancombc2(
  data = data_phyloseq,
  fix_formula = 'substrate + timepoint',
  p_adj_method = 'holm',
  group = 'timepoint',
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  verbose = TRUE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
)


# ===== Check ANCOMBC2 pairwise results

res_pair <- output$res_pair

# Uncomment code below if you want to save pairwise results table
# write.table(
#   res_pair, 
#   file='tables_and_figures/asv/ancom/timepoint/ancom_results_timepoint.tsv', 
#   sep='\t', 
#   row.names = FALSE
# )

# The long line of code below parses the pairwise results table to view
#   only the features that are differentially abundant in at least
#   1 pair of sample group
# In the pairwise results table, you may see some columns showing two
#   sample groups (e.g. timepoint48h_timepoint24h). This means that
#   the statistics (e.g. log fold change) for the 48h sample group for
#   a certain feature is calculated with respect to 24h.
#   For cases where there is no 'baseline' group (e.g. timepoint24h),
#   that  means that it is compared to the first category in the
#   metadata variable of interest. In this case, timepoint24h
#   implicitly means timepoint24h versus timepoint12h
# Modify the code below as needed. You can also refer to the tutorial
#   linked above for a more detailed step-by-step discussion.
# Briefly, the code below does the following:
#   (1) filter - select only the features that are differentially abundant
#       in at least 1 pairing
#   (2) mutate - round off log fold change (lfc) and p-value
df_fig_pair1 <- res_pair %>%
  filter(
    diff_timepoint24h == 1 | 
    diff_timepoint48h == 1 | 
    diff_timepoint7d == 1 |
    diff_timepoint48h_timepoint24h == 1 |
    diff_timepoint7d_timepoint24h == 1 |
    diff_timepoint7d_timepoint48h == 1
  ) %>%
  mutate(
    lfc1 = ifelse(
      diff_timepoint24h == 1, 
      round(lfc_timepoint24h, 2), 0
    ),
    lfc2 = ifelse(
      diff_timepoint48h == 1, 
      round(lfc_timepoint48h, 2), 0
    ),
    lfc3 = ifelse(
      diff_timepoint7d == 1, 
      round(lfc_timepoint7d, 2), 0
    ),
    lfc4 = ifelse(
      diff_timepoint48h_timepoint24h == 1, 
      round(lfc_timepoint48h_timepoint24h, 2), 0
    ),
    lfc5 = ifelse(
      diff_timepoint7d_timepoint24h == 1, 
      round(lfc_timepoint7d_timepoint24h, 2), 0
    ),
    lfc6 = ifelse(
      diff_timepoint7d_timepoint48h == 1, 
      round(lfc_timepoint7d_timepoint48h, 2), 0
    ),
    stat1 = ifelse(
      diff_timepoint24h == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint24h, 2), q_timepoint24h), NA
    ),
    stat2 = ifelse(
      diff_timepoint48h == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint48h, 2), q_timepoint48h), NA
    ),
    stat3 = ifelse(
      diff_timepoint7d == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint7d, 2), q_timepoint7d), NA
    ),
    stat4 = ifelse(
      diff_timepoint48h_timepoint24h == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint48h_timepoint24h, 2), q_timepoint48h_timepoint24h), NA
    ),
    stat5 = ifelse(
      diff_timepoint7d_timepoint24h == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint7d_timepoint24h, 2), q_timepoint7d_timepoint24h), NA
    ),
    stat6 = ifelse(
      diff_timepoint7d_timepoint48h == 1, 
      sprintf('%.2f (%.2e)', round(lfc_timepoint7d_timepoint48h, 2), q_timepoint7d_timepoint48h), NA
    )
  )

df_fig_pair1 # Print results of differentially abundant features

# Melt sub-table of log fold changes
lfc_pivot <- df_fig_pair1 %>%
  pivot_longer(
    cols = lfc1:lfc6,
    names_to = 'group',
    values_to = 'value'
  ) %>%
  arrange(taxon)

# Subset columns
lfc_pivot <- lfc_pivot[c('taxon', 'group', 'value')]

# Melt sub-table of p-values
stat_pivot <- df_fig_pair1 %>%
  pivot_longer(
    cols = stat1:stat6,
    names_to = 'group',
    values_to = 'stat'
  )

# Subset columns
stat_pivot <- stat_pivot[c('taxon', 'group', 'stat')]

# Rename lfc categories
# Check the long mutate function in the long code above to know
#   which lfc category corresponds to which group comparison
lfc_pivot$group <- recode(
  lfc_pivot$group,
  'lfc1' = '24h - 12h',
  'lfc2' = '48h - 12h',
  'lfc3' = '7d - 12h',
  'lfc4' = '48h - 24h',
  'lfc5' = '7d - 24h',
  'lfc6' = '7d - 48h'
)

# Rename stat categories
# Check the long mutate function in the long code above to know
#   which lfc category corresponds to which group comparison
stat_pivot$group <- recode(
  stat_pivot$group,
  'stat1' = '24h - 12h',
  'stat2' = '48h - 12h',
  'stat3' = '7d - 12h',
  'stat4' = '48h - 24h',
  'stat5' = '7d - 24h',
  'stat6' = '7d - 48h'
)

# Unify lfc_pivot and stat_pivot into a single table
summary_df <- merge(lfc_pivot, stat_pivot, by=c('taxon', 'group'))


# ===== Rename ASV IDs to Taxonomy names =====

get_rightmost_non_na <- function(x) {
  non_na_indices <- which(!is.na(x))
  if (length(non_na_indices) > 0) {
    return(x[tail(non_na_indices, 1)])
  } else {
    return(NA)
  }
}

taxon_name <- apply(taxonomy_table[summary_df$taxon, ], 1, get_rightmost_non_na)
taxon_name

summary_df$taxon_name <- revalue(summary_df$taxon, taxon_name)


# ===== Plot heatmap =====

# Define the LFC lowest value
lo <- ifelse(
  abs(floor(min(summary_df$value))) >= ceiling(max(summary_df$value)),
  floor(min(summary_df$value)),
  -ceiling(max(summary_df$value))
)

# Define the highest LFC value
up <- ifelse(
  abs(floor(min(summary_df$value))) >= ceiling(max(summary_df$value)),
  -floor(min(summary_df$value)),
  ceiling(max(summary_df$value))
)

# Generate plot
# Uncomment geom_text if you want p-values to be written in the graph as well
fig_pair <- ggplot(summary_df, aes(x=group, y=taxon, fill=value)) +
  geom_tile(color='black') +
  scale_fill_gradient2(
    low='blue',
    high='red',
    mid='white',
    na.value='white',
    midpoint=0,
    limit=c(lo, up),
    name=NULL
  ) +
  scale_color_identity(guide='none') +
  labs(
    x=NULL, 
    y=NULL, 
    title='Log fold changes'
  ) +
  scale_y_discrete(labels=taxon_name) +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5))
  # geom_text(
  #   aes(
  #     group,
  #     taxon,
  #     label=stat,
  #   ),
  #   size=2
  # )

# Show plot
fig_pair

# Uncomment code below if you want to save the plot
# Adjust parameters height and width
# ggsave(
#   'tables_and_figures/asv/ancom/timepoint/ASV_LFC_as_taxon_names_corrected.svg',
#   fig_pair,
#   width = 2300,
#   height = 5000,
#   units = 'px',
#   bg = 'white'
# )


