=======================================
DESCRIPTION
=======================================
-This folder is a storage for the raw and processed data that can be used for the demonstration of the workflows in the jupyter notebooks inside the metabarcoding_resources folder

=======================================
CONTENT
=======================================
(1) asv_table_and_sequences - feature (ASV) table and sequences, and taxonomic assignments produced by the commands outlined in the denoising pipeline; used as inputs to the Amplicon_Feature_Filtering.ipynb notebook for demo

(2) otu_table_and_sequences - feature (OTU) table and sequences, and taxonomic assignments produced by the commands outlined in the clustering pipeline; used as inputs to the Amplicon_Statistical_Analyses.ipynb notebook for demo

(3) output_data/asv_denoising/ - all data produced by the denoising pipeline; if you are having trouble in running a certain step, you can simply copy and paste the output of that step from this folder so that you could proceed to the next step

(4) output_data/otu_clustering/ - all data produced by the clustering pipeline; if you are having trouble in running a certain step, you can simply copy and paste the output of that step from this folder so that you could proceed to the next step

(5) picrust2_input_data/ - feature table and representative sequences to be used in the PICRUSt2_Pipeline.ipynb notebook. Data was taken from Raes et al. (2021): 

https://zenodo.org/record/4567694#.YhnG6ehBxPZ

(6) output_data/picrust2/ - all data produced by the picrust2 pipeline; if you are having trouble in running a certain step, you can simply copy and paste the output of that step from this folder so that you could proceed to the next step

(7) metadata.txt - metadata file used in the denoising and clustering workflows

(8) metadata-day-2.txt - metadata file used in the Amplicon_Statistical_Analyses.ipynb notebook

(9) data-links.txt - ftp links to the input data for Amplicon_Clustering_Pipeline.ipynb and Amplicon_Denoising_Pipeline.ipynb