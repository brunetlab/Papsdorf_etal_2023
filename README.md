# All code associated with Papsdorf et. al. 2023
## Scripts:
Scripts 1-8 are code analyzing the lipidomics data and generate the corresponding figures and tables. The lipidomics data contains conditions, that are not used in the current manuscript (PA and PRX). These are omitted from the plots.

Script 9 generates the GO-term heatmap

Script 10 and 11 calculate the lipid droplet peroxisome screen statistics

## Raw data:
Contains all the raw data that is necessary to run the scripts

## Table of scripts and corresponding figures
If a script relies on output from previous scripts, this is indicated in the table below. The corresponding figure panels are listed below.

|SCRIPT NAME|REQUIRED OUTPUTS FROM <br>OTHER SCRIPTS TO RUN| OUTPUT|FIGURE PANELS/ TABLES|
  ---|---|---|---
  |1_LIPIDS_PREPROCESSING_IS_<br>NORMALIZATION_20221010|no|1_Lipids_Preprocessed_Normalized.csv||
  |||1_Lipids_Preprocessed_Normalized.Rdata||
  |2_LIPIDS_IMPUTATION_CV_FILTERING|1_Lipids_Preprocessed_<br>Normalized.Rdata|2_Lipids_Imputed_CVfiltered.csv |Extended Data Table 3|
  |||2_Lipids_Imputed_CVfiltered.Rdata||
  |3_LIPIDS_PCA_20221010|2_Lipids_Imputed_CVfiltered.Rdata|3_PCA.pdf |Figure 4b|
  |4_Lipid_Classes_Triglycerides_20221010|	2_Lipids_Imputed_CVfiltered.Rdata	|4_Classes.csv|Data Table 1|
  |||	4_Class_TG.pdf|Extdended Data Figure 4a|
  |5_SFA_MUFA_PUFS_All_Lipids	|2_Lipids_Imputed_CVfiltered.Rdata|5_SFA_MUFA_PUFA_All_Lipids.csv| Data Table 1|	
  |||5_SFA-MUFA-PUFA_All_Lipids.pdf|Figure 4c|
  |6_SFA_MUFA_PUFA_Phospholipids	|2_Lipids_Imputed_CVfiltered.Rdata|6_SFA_MUFA_PUFA_Phospholipids.csv|	Data Table 1 |
  |||6_MUFAtoPUFA_Phospholipids.pdf|Figure 4d|
  |7_Ether_Lipid_Abundance_20221010|	2_Lipids_Imputed_CVfiltered.Rdata|	7_Ether_Lipids.csv	| Data Table 1 |
  |||7_Ether_Lipids.pdf|Figure 4e|
  |8_Peroxidation_Index|	2_Lipids_Imputed_CVfiltered.Rdata	|8_Peroxidation_Index_All_Lipids.csv|	Data Table 1|
  |||8_PI_All_Lipids.pdf|Figure 4g|
  |9_Shared_GO_terms_upregulated_Greer_Han|	no	|9_Shared_GO_terms.csv	|Data Table 1|
  |||9_Shared_Go_term_heatmap|Figure 5b|
  |10_LD_stats_20221201|	no	|10_Lipid_droplet_stats_20221201.csv	|Data Table 1, coloring of significance in Extended Data Figure 5n|
  |11_P_stats_20221201|	no	|11_Peroxisome_stats_20221201.csv	|Data Table 1, coloring of significance in Extended Data Figure 5o|
  
