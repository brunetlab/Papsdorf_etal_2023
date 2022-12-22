######## Preprocessing and Internal Standard Normalization #################################################

# General info:
#               For the manuscript Papsdorf_2022 only the conditions control (EV)and ash-2 RNAi (ASH) are plotted
#               Analysis is performed on all four conditions (EV, ASH, PRX-5 and PRX-5/ASH2 RNAi) prepared in sextuplicates. 
#               Internal standard: 2uL equiSPLASH

# Data: Raw data 
# Approach: 1)	Load data, select the correct ionization mode
          # 2)	Filter data lipids with 0 signal and aggregate signals from the same lipid
          # 3)	Load internal standards, select the correct ionization mode
          # 4)  Normalize using internal standard (equiSPLASH 2uL-> for molar values see table Standard concentration in Raw data).
          #     CL, DLCL, MLCL and PA are quantified using PG standard, LPI, LPS are quantified using LPE standard, and MG using DG standard due to molecular similarity.
          # 5)	Subtract blank and discard rows with more than 50% NA
          # 6)  Normalize: dividing by protein concentration (nmol lipid/mg protein)
          # 7)  Normalize: dividing by sample median * global median


#Output: 1_Lipids_Preprocessed_Normalized.Rdata  (in ng/mg protein)

#Packages:
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location

#Create directories for figure creation

dir.create("../Output_Data")
dir.create("../Output_Figures")

rm(list=ls())


################## 1)	Load data, select the correct ionization mode #################################################

#Load negative mode data
table_neg<-read.csv("../Raw_data/210223_KP_neg.csv", fill = TRUE, header = TRUE, row.names=NULL, stringsAsFactors = FALSE)[,1:38]
table_neg<-table_neg[,-c(1,3,8:12,37)] #Delete columns that are not needed for analysis (LipidIon/fatty acids merged/CalcMZ to GroupTopPos.c./QC MSMS)

#Select the correct ionization mode for per lipid. Negative mode data contains "CL","DLCL","LPC","LPE","LPI","LPS","MLCL","PA","PC","PE","PG","PI","PS")
table_neg<-table_neg %>%
  filter(Class %in% c("CL","DLCL","LPC","LPE","LPI","LPS","MLCL","PA","PC","PE","PG","PI","PS"))

#Load positive mode data
table_pos<-read.csv("../Raw_data/210223_KP_pos.csv", fill = TRUE, header = TRUE, row.names=NULL, stringsAsFactors = FALSE)[,1:38]
table_pos<-table_pos[,-c(1,3,8:12,37)] #Delete columns that are not needed for analysis (LipidIon/fatty acids merged/CalcMZ to GroupTopPos.c./QC MSMS)

#Select the correct ionization mode for per lipid. Positive mode data contains "Cer","ChE","DG","MG","SM","TG"
table_pos<-table_pos %>%
  filter(Class %in% c("Cer","ChE","DG","MG","SM","TG"))

# Bind tables
table<-bind_rows(table_neg,table_pos)


################## 2)	Filter lipids with 0 signal and aggregate signals from the same lipid #############################

table <- table %>%
  filter_at(vars(6:29), any_vars((.) != 0))%>%  # Delete lipids with 0 signal
  group_by(Class,FA1,FA2,FA3,FA4) %>% 
  summarise_all(sum,na.rm = FALSE)  # Aggregate


################## 3)	Load internal standards, select the correct ionization mode ######################################

#load negative mode internal standard
is_neg<-read.csv("../Raw_data/210223_KP_neg.csv", fill = TRUE, header = FALSE, row.names=NULL, stringsAsFactors = FALSE)[1:29,40:48]

is_neg[1,]<-gsub("-","",is_neg[1,]) #delete the - in the column names

is_neg[1,]<-gsub("d","",is_neg[1,]) #delete the d in the column names

is_neg[,1]<-gsub("-",".",is_neg[,1]) #delete the - in the row names

colnames(is_neg)<- is_neg[1,] #set colnames

is_neg<-is_neg[-c(1:3,5),] #delete the QCs

colnames(is_neg)[1]<-"Sample"

#load positive mode internal standard
is_pos<-read.csv("../Raw_data/210223_KP_pos.csv", fill = TRUE, header = FALSE, row.names=NULL, stringsAsFactors = FALSE)[1:28,40:45]

is_pos[1,]<-gsub("d","",is_pos[1,]) #delete the - in the column names

is_pos[,1]<-gsub("-",".",is_pos[,1]) #delete the d in the column names

colnames(is_pos)<- is_pos[1,] #set colnames

is_pos<-is_pos[-c(1:3),] #delete the QCs

colnames(is_pos)[1]<-"Sample"

#make internal standards table
is_table<-bind_cols(is_neg[,-7],is_pos[,-1]) # Delete names of positive standards and Ceramide standard (Cer) form negative mode.

# Add standards for quantification by using molecular similarity
is_table<-is_table %>%
  mutate(CL = PG) %>% #annotate the same values from CL to PG etc
  mutate(DLCL = PG) %>%
  mutate(MLCL = PG) %>%
  mutate(LPI = LPE) %>%
  mutate(LPS = LPE) %>%
  mutate(MG = DG)


################## 4) Normalize using IS (2ul equiSPLASH 2ul -> load IS Stanfard concentrations in nmol for normalization ) ######################################

is_molarity<-read.csv("../Raw_data/IS Standard concentrations nmol in 2 ul.csv", fill = TRUE, header = FALSE, row.names=NULL, stringsAsFactors = FALSE)
row.names (is_molarity)<-is_molarity[,1]
is_molarity<-is_molarity [-c(1)]

table_norm<-table%>%
  dplyr::select(c(1:5)) # make new table with lipid names

for (j in sort(is_table$Sample)) { 
  
  temp_IS<-is_table%>%filter(Sample == j) 
  temp_tab<-table%>%dplyr::select(c(1:5,j)) 
  class_norm<-tibble() 
  
  
  for (i in colnames(is_table)[2:ncol(is_table)])  { 
    IS<-temp_IS%>%dplyr::select(i) 
    temp_class<-temp_tab%>%filter(Class == i) 
    norm_factor<-is_molarity %>% filter(row.names(is_molarity) %in% c(i)) #filter single value based on IS from is_molarity
    y=norm_factor[1,1]
    temp_class[,ncol(temp_class)]<-temp_class[,ncol(temp_class)]/as.numeric(IS)*y #in nmol. Normalization step by dividing lipid value/internal standard
    class_norm<-bind_rows(class_norm,temp_class) #adds the standardized values to normalized table
  }
   table_norm<-full_join(table_norm,class_norm) 
}

rm("i")
rm("j")
rm("y")


################## 5) Subtract blank ###########################################################################

table_norm[,6:ncol(table_norm)]<-table_norm[,6:ncol(table_norm)]-((table_norm$blank)*3) # Everything that is lower than 3 times the blank is discarded
table_norm<-table_norm %>% #delete the blank from the dataset
  select(-blank)

# Remove rows with more than 50% <=0 values in the dataset
table_norm[,6:29][table_norm [,6:29] <= 0] <- NA # <0 =  NA
table_norm<-table_norm[which(rowMeans(!is.na(table_norm [,6:29])) > 0.5), ] # Remove rows with more than 50% NA


################## 6)	Normalize by protein concentration #######################################################

# Load protein concentration and normalize
prot_cc_tab<-read.csv("../Raw_data/210223_KP_protein_concentration_mgml.csv", header = TRUE) 

prot_cc_tab<- prot_cc_tab %>% 
  select(sort(names(prot_cc_tab)))

prot_cc_tab<-prot_cc_tab*50/1000 # (mg of protein in sample) 

table_norm_p<-as.data.frame(table_norm)

table_norm_p[,6:ncol(table_norm_p)]<-table_norm_p[,6:ncol(table_norm_p)]/as.list(as.numeric(prot_cc_tab[1,])) # (nmol of lipid/mg of protein) 


# Eliminate brackets
for (i in 2:5) {
  table_norm_p[,i] <- gsub("\\(|\\)","",table_norm_p[,i])
  table_norm_p[,i]  <- gsub("\\)","",table_norm_p[,i])
}


################## 7)	Normalize by lipid median #############################################################

# Calculate normalization factor = median of each sample * global median
table_norm_med<-table_norm_p
median_table <- apply(table_norm_med[,6:ncol(table_norm_med)], 2, median, na.rm=T) #make median for each column (except names)
global_median = median(median_table)


# Normalize to median
table_norm_med[,6:ncol(table_norm_med)] <- sweep(table_norm_med[,6:ncol(table_norm_med)], 2, median_table, "/")
table_norm_med[,6:ncol(table_norm_med)] <- sweep(table_norm_med[,6:ncol(table_norm_med)], 2, global_median, "*") #output in nmol/mg protein


# Calculate max sum for all lipids to use for normalization purposes for lipid oxidation
table_max_sum <- table_norm_p
table_max_sum<-as.data.frame(apply(table_max_sum[,6:ncol(table_max_sum)], 2, sum, na.rm=T)) #make median for each column (except names)
write.csv(table_max_sum, "../Output_Data/1_Lipids_Preprocessed_Max_SUM.csv")

################## Save

save(table_norm_med, file = paste0("../Output_Data/1_Lipids_Preprocessed_Normalized", ".Rdata")) #output in nmol/mg protein
write.csv(table_norm_med, "../Output_Data/1_Lipids_Preprocessed_Normalized.csv")



# sessionInfo(package = NULL)
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.5     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3     tibble_3.1.1    ggplot2_3.3.3   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8       cellranger_1.1.0 pillar_1.6.0     compiler_3.6.3   dbplyr_2.1.1     tools_3.6.3      jsonlite_1.7.2   lubridate_1.7.10 lifecycle_1.0.0 
# [10] gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.10     reprex_2.0.0     cli_2.5.0        DBI_1.1.0        rstudioapi_0.13  haven_2.3.1      xml2_1.3.2      
# [19] withr_2.4.2      httr_1.4.2       fs_1.5.0         generics_0.1.0   vctrs_0.3.7      hms_1.0.0        grid_3.6.3       tidyselect_1.1.0 glue_1.4.2      
# [28] R6_2.5.0         fansi_0.4.1      readxl_1.3.1     modelr_0.1.8     magrittr_2.0.1   backports_1.2.0  scales_1.1.1     ellipsis_0.3.2   rvest_1.0.0     
# [37] assertthat_0.2.1 colorspace_1.4-1 utf8_1.1.4       stringi_1.5.3    munsell_0.5.0    broom_0.7.6      crayon_1.4.1
