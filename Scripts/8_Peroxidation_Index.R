######## Calculate Peroxidation Index #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#      normalized, imputed, CV_filtered
#      Only the samples EV (EV) and ash-2 RNAi(ASH) are plotted

# Approach: 1) Fatty acyl chain abundance in all lipids (individual fatty acyl concentration matches the corresponding lipid). Odd chain fatty acids and fatty acids with d, O and t modifications are not analyzed.
#           2) Categorize by number of insaturations
#           3) Calculate Peroxidation index with the following formula PIn=0.025×(percentage of monoenoics)+1×(percentage of dienoics)+2×(percentage of trienoics)+4×(percentage of tetraenoics)+6×(percentage of pentaenoics)+8×(percentage of hexaenoics) 
#               Formula from here : https://journals.biologists.com/jeb/article/223/11/jeb224063/223591/Membrane-peroxidation-index-and-maximum-lifespan
#           4) Fold change and stats 
#           5) Plot 

# Goal: Plot Peroxidation index

# Output:  Pdf: 17_Peroxidation_index_all_lipids.pdf
#          Csv: 17_Peroxidation_index_all_lipids.csv   



rm(list=ls())
library(ggplot2)
library(tidyverse)
library(ggpubr)

options(stringsAsFactors = FALSE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location


load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered



######## 1) Fatty acyl abundance in all lipids #########################################################################################

########FA count. List with individual fatty acyl side chains listed with the corresponding lipid concentration in one table

FA_table<-final_table[,c(2,6:ncol(final_table))] #make a table with all FA1

FA_temp<-final_table[,c(3,6:ncol(final_table))] ##make a table with all FA2

colnames(FA_temp)[1]<-c("FA1") #adjust colnames so that they match for bind_rows

FA_table<-bind_rows(FA_table,FA_temp) #merge

FA_temp<-final_table[,c(4,6:ncol(final_table))] #repeat step above with FA3 

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_temp<-final_table[,c(5,6:ncol(final_table))] #repeat step above with FA4

colnames(FA_temp)[1]<-c("FA1")

FA_table<-bind_rows(FA_table,FA_temp)

FA_table_1<-bind_rows(FA_table,FA_temp)

colnames(FA_table_1)[1]<-"FA"


#Filter out odd chain fatty acids, O,t,d modifications 

FA_table_1<- filter(FA_table_1, !grepl("t",FA))

FA_table_1<- filter(FA_table_1, !grepl("d",FA))

FA_table_1<- filter(FA_table_1, !grepl("O",FA))

uns_FA_table<- FA_table_1 %>% 
  select(1:25) %>%
  filter(!grepl("6:0",FA)) %>%
  filter(!grepl("1:",FA)) %>%
  filter(!grepl("3:",FA)) %>%
  filter(!grepl("5:",FA)) %>%
  filter(!grepl("7:",FA)) %>%
  filter(!grepl("9:",FA))

######## 2) Categorize fatty acids by number of insaturations ##########################################################################

#Table with the number of unsaturations (ether-linked fatty acids are kept)

uns_FA_table<-uns_FA_table %>%
  mutate(FA = str_replace_all(FA, "e", "")) %>% #replace all the "e"s with spaces
  mutate(Uns = sub(".*:", "", FA)) %>% #delete the info before the ":"
  select(Uns,2:25)


#Aggregate

uns_FA_table<-aggregate(uns_FA_table[,-1], by = list(uns_FA_table$Uns), FUN = sum) #summarize intensity of the same fatty acid unsaturations

uns_FA_table<-uns_FA_table[-1,]#delete the first row with the empty FA loadings 

uns_FA_table <-as.data.frame(lapply(uns_FA_table,as.numeric)) #make numeric


########3) Calculate peroxidation index #################################################

#Calculate % of insaturations

uns_FA_table1 <- uns_FA_table %>% mutate_all(sum) #summarize to calculate percentage

uns_FA_table2<-100/uns_FA_table1*uns_FA_table #calculate %

uns_FA_table2[,1]<-uns_FA_table[,1] # change to correct name (number of insaturations)

uns_FA_table3<-uns_FA_table2[-1,] #remove row with 0 double bonds (not used for calculations)


#Calculate Peroxidation Index

uns_FA_table3$x <- c(0.025, 1, 2, 4, 6, 8) #make column with factor for multiplication (x). Double check what is the maximum of double bonds and adjust

uns_FA_table4<-as.data.frame(lapply(uns_FA_table3, `*`, uns_FA_table3$x)) #multiply factor with x

uns_FA_table5 <- uns_FA_table4 %>% mutate_all(sum) #summarize to calculate percentage#sum

uns_FA_table5 <-uns_FA_table5[1,]


######## 4) Fold change and stats ###########################################################################################################

SMP_table<-uns_FA_table5[,-1]
SMP_table<-SMP_table[,-25]

#Mean, Coefficient of Variation and Fold Change

SMP_table<-SMP_table %>% 
  mutate(ASH = apply(SMP_table %>% select(1:6), 1, mean ))%>% 
  mutate(EV = apply(SMP_table %>% select(7:12), 1, mean ))%>%
  mutate(PA = apply(SMP_table %>% select(13:18), 1, mean ))%>% 
  mutate(PRX = apply(SMP_table %>% select(19:24), 1, mean ))%>%
  mutate(cvASH = apply(SMP_table %>% select(1:6), 1, sd )/ASH)%>% 
  mutate(cvEV = apply(SMP_table %>% select(7:12), 1, sd )/EV)%>%
  mutate(cvPA = apply(SMP_table %>% select(13:18), 1, sd )/PA)%>% 
  mutate(cvPRX = apply(SMP_table %>% select(19:24), 1, sd )/PRX)%>%
  mutate(FC_EVvsASH = (ASH/EV))%>%
  mutate(FC_PRXvsPA = (PRX/PA))

#Wilcoxon

for (i in 1:nrow(SMP_table)) {
  
  SMP_table$pval_EVvsASH[i]<-wilcox.test(as.numeric(SMP_table[i,1:6]),as.numeric(SMP_table[i,7:12]))$p.value #to compare ASH vs EV
  SMP_table$pval_PRXvsPA[i]<-wilcox.test(as.numeric(SMP_table[i,19:24]),as.numeric(SMP_table[i,13:18]))$p.value #to compare PRX-ASH with PRX
  SMP_table$pval_EVvsPA[i]<-wilcox.test(as.numeric(SMP_table[i,7:12]),as.numeric(SMP_table[i,13:18]))$p.value #to compare EV with PRX-ASH
  SMP_table$pval_EVvsPRX[i]<-wilcox.test(as.numeric(SMP_table[i,7:12]),as.numeric(SMP_table[i,19:24]))$p.value #to compare EV with PRX
}

#Adjust p-val

SMP_table<- SMP_table %>%
  mutate(FDR_EVvsASH = p.adjust(pval_EVvsASH,method = "BH"))%>%
  mutate(FDR_PRXvsPA = p.adjust(pval_PRXvsPA,method = "BH")) %>%
  mutate(FDR_EVvsPA = p.adjust(pval_EVvsPA,method = "BH"))%>%
  mutate(FDR_EVvsPRX = p.adjust(pval_EVvsPRX,method = "BH"))


#Add names

SMP_table$Peroxidation_index = "Peroxidation_index" 
SMP_table <- SMP_table %>%
  select(Peroxidation_index, everything()) #move SMP to the front
SMP_table$Peroxidation_index <-rownames(SMP_table)


#Save
write.csv(SMP_table, '../Output_Data/8_Peroxidation_Index_All_Lipids.csv')


########Organize data for plot

tSMP_table<-as.data.frame(t(SMP_table[,2:13])) %>%
  mutate(Group = NA)

tSMP_table[c(1:6),ncol(tSMP_table)]<- "ASH"

tSMP_table[c(7:12),ncol(tSMP_table)]<- "EV"

tSMP_table$Group <- factor(tSMP_table$Group , levels=c("EV", "ASH")) #order the groups

ptSMP_table<-pivot_longer(tSMP_table, cols = c(1:(ncol(tSMP_table)-1)), values_to = "Val")



######## 5) Plot         ###################################################################################################################

######## Plot Peroxidation index        

plot<- ptSMP_table

PI<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0, outlier.size=1.5)+
  geom_point(data=plot, aes(x=Group, y=Val), position=position_jitter(width=0.15)) +
  xlab("")+
  ylab("Peroxidation index all lipids")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5), plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values = c("EV" = "grey55",
                               "ASH" = "red1"))

plot(PI)

########Save Plot

ggsave("../Output_Figures/8_PI_All_Lipids.pdf", PI, width=0.7, height=1.4, units="in", scale=3,useDingbats=FALSE) 

# 
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
#   [1] ggpubr_0.4.0.999  data.table_1.13.6 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.5       purrr_0.3.4       readr_1.4.0       tidyr_1.1.3      
# [9] tibble_3.1.1      ggplot2_3.3.3     tidyverse_1.3.1  
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
