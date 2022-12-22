######## Screen stats Peroxisomes #################################################

# Goal: Stats for lipid/peroxo screen 

# Data: Normalized lipid droplets and peroxisome counts

# Approach: (1) Calculate exact p-values and p-adjust

# Output:   Csv: 11_P_stats_20221201.csv 


#Packages
library(pheatmap)
library(tidyverse)
library (ggplot2)

rm(list=ls())


# Set wd to current file location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

dfP<- read.csv ('../Raw_data/Normalized P numbers 20220411.csv', stringsAsFactors = FALSE)
dfP<- as.data.frame(dfP)

###################################################
##########       Peroxisomes            ###########
###################################################



################ Calculate exact p-values peroxisomes  ################################################################################

x=as.numeric (dfP[,"Control"]) #to compare each condition to the control (x)

#calculate p-values in loop

p_val_P<- NULL

for (i in 1:ncol(dfP))
{
  p_value<-wilcox.test(x,dfP[,i])[["p.value"]] 
  {
    p_val_P[[length(p_val_P) + 1]] <- p_value
  }
}

p_val_P<-as.data.frame(p_val_P)
rownames(p_val_P)<-colnames(dfP)

colnames(p_val_P)<-c("p_val_P")

# check for 2 values-->correct
# 
# i=5 #(2=ash-2 and 5=plin-1)
# p_value<-(wilcox.test(x,dfP[,i])[["p.value"]])
# p_value #ash-2 = 1.359617e-30 , 0.00762936
#correct-results overlap with table
# 
# rm(i)


################ Calculate padjust P

p_adj_P<- NULL
#p.adjust(class_ntable$pval,method = "BH")

for (i in 1:ncol(p_val_P))
{
  p_adj<-p.adjust(p_val_P$p_val_P,method = "BH")
  {
    p_adj_P[[length(p_adj_P) + 1]] <- p_adj
  }
}

p_adj_P<-as.data.frame(p_adj_P)
row.names(p_adj_P)<-colnames(dfP)
colnames(p_adj_P)<-"p_adj_P"


################ Merge and sort alphabetically except for first three rows
P_stats <-cbind(p_val_P, p_adj_P)

Psub <- P_stats[4:nrow(P_stats),]
Psub<-Psub[order(row.names(Psub)), ]

P_stats_final <- rbind(P_stats[1:3,], Psub)

write.csv(P_stats_final, "../Output_Data/11_Peroxisome_stats_20221201.csv")

# > sessionInfo(package = NULL)
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
#   [1] RColorBrewer_1.1-3 pheatmap_1.0.12    ggpubr_0.4.0.999   data.table_1.13.6  forcats_0.5.1      stringr_1.4.0      dplyr_1.0.5       
# [8] purrr_0.3.4        readr_1.4.0        tidyr_1.1.3        tibble_3.1.1       ggplot2_3.3.3      tidyverse_1.3.1   
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.0 haven_2.3.1      carData_3.0-4    colorspace_1.4-1 vctrs_0.3.7      generics_0.1.0   utf8_1.1.4       rlang_0.4.10    
# [9] pillar_1.6.0     foreign_0.8-76   glue_1.4.2       withr_2.4.2      DBI_1.1.0        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
# [17] lifecycle_1.0.0  munsell_0.5.0    ggsignif_0.6.1   gtable_0.3.0     cellranger_1.1.0 zip_2.1.1        rvest_1.0.0      rio_0.5.26      
# [25] labeling_0.4.2   curl_4.3         fansi_0.4.1      broom_0.7.6      Rcpp_1.0.8       scales_1.1.1     backports_1.2.0  jsonlite_1.7.2  
# [33] abind_1.4-5      farver_2.0.3     fs_1.5.0         hms_1.0.0        digest_0.6.27    openxlsx_4.2.3   stringi_1.5.3    rstatix_0.7.0   
# [41] cowplot_1.1.1    grid_3.6.3       cli_2.5.0        tools_3.6.3      magrittr_2.0.1   car_3.0-10       crayon_1.4.1     pkgconfig_2.0.3 
# [49] ellipsis_0.3.2   MASS_7.3-53      xml2_1.3.2       reprex_2.0.0     lubridate_1.7.10 assertthat_0.2.1 httr_1.4.2       rstudioapi_0.13 
# [57] R6_2.5.0         compiler_3.6.3  



