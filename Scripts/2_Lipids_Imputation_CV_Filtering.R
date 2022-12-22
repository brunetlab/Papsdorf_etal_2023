######## Imputation and CV Filtering #################################################

# Data: 1_Lipids_Preprocessed_Normalized (Output_Data)

# Approach: 1) Impute 330 NAs by assigning a random value based on the bottom 5% for the corresponding lipid. 
          # 2) Filtering for a coefficient of variance <0.5

         
#Output: 2_Lipids_Imputed_CVfiltered.Rdata      (in (nmol/mg protein))
#        2_Lipids_Imputed_CVfiltered.csv

#Packages:
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location
rm(list=ls())

load(file = "../Output_Data/1_Lipids_Preprocessed_Normalized.Rdata")
df.impute <-table_norm_med 


################ 1) Imputation for NA values ####################################################################################

length(df.impute[,c(1,6:ncol(df.impute))][is.na(df.impute[,c(1,6:ncol(df.impute))])]) #330 NAs in dataset

#Impute
for (i in 1:nrow(df.impute[,c(1,6:ncol(df.impute))]))
     {
  b <- df.impute[,c(1,6:ncol(df.impute))][i, which(is.na(df.impute[,c(1,6:ncol(df.impute))][i,])==T)]
  df.impute[,c(1,6:ncol(df.impute))][i, which(is.na(df.impute[,c(1,6:ncol(df.impute))][i,])==T)]<- rnorm(n = length(b),
                                                               mean = quantile(as.numeric(df.impute[,c(1,6:ncol(df.impute))][i,]),na.rm = T, 0.05),
                                                               sd = 0.3*quantile(as.numeric(df.impute[,c(1,6:ncol(df.impute))][i,]),na.rm = T, 0.05))
}
length(df.impute[,c(1,6:ncol(df.impute))][is.na(df.impute[,c(1,6:ncol(df.impute))])]) #0

#Check
b <- setdiff(table_norm_med, df.impute)
c <- setdiff(df.impute, table_norm_med)
length(which(!b==c)) ##0 false value


################ 2) CV calculation and Filtering ####################################################################################

df.impute_2<- df.impute %>%
  mutate(CV_ASH = apply(df.impute %>% select(6:11), 1, sd )/apply(df.impute %>% select(6:11), 1, mean )) %>%
  mutate(CV_EV = apply(df.impute %>% select(12:17), 1, sd )/apply(df.impute %>% select(12:17), 1, mean ))%>%
  mutate(CV_PRX = apply(df.impute %>% select(24:29), 1, sd )/apply(df.impute %>% select(24:29), 1, mean )) %>%
  mutate(CV_PRX_ASH = apply(df.impute %>% select(18:23), 1, sd )/apply(df.impute %>% select(18:23), 1, mean ))


df.impute_CVfiltered<- df.impute_2 %>% 
  filter(CV_ASH < 0.5) %>%
  filter(CV_EV < 0.5) %>%
  filter(CV_PRX < 0.5) %>%
  filter(CV_PRX_ASH < 0.5) %>%
  select(-CV_EV,-CV_ASH,-CV_PRX_ASH,-CV_PRX)


################## Save

save(df.impute_CVfiltered, file = paste0("../Output_Data/2_Lipids_Imputed_CVfiltered", ".Rdata"))
write.csv(df.impute_CVfiltered, "../Output_Data/2_Lipids_Imputed_CVfiltered.csv")


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
#   [1] Rcpp_1.0.8       cellranger_1.1.0 pillar_1.6.0     compiler_3.6.3   dbplyr_2.1.1     tools_3.6.3      lubridate_1.7.10 jsonlite_1.7.2   lifecycle_1.0.0 
# [10] gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.10     reprex_2.0.0     cli_2.5.0        DBI_1.1.0        rstudioapi_0.13  haven_2.3.1      xml2_1.3.2      
# [19] withr_2.4.2      httr_1.4.2       fs_1.5.0         generics_0.1.0   vctrs_0.3.7      hms_1.0.0        grid_3.6.3       tidyselect_1.1.0 glue_1.4.2      
# [28] R6_2.5.0         fansi_0.4.1      readxl_1.3.1     modelr_0.1.8     magrittr_2.0.1   backports_1.2.0  scales_1.1.1     ellipsis_0.3.2   rvest_1.0.0     
# [37] assertthat_0.2.1 colorspace_1.4-1 utf8_1.1.4       stringi_1.5.3    munsell_0.5.0    broom_0.7.6      crayon_1.4.1
