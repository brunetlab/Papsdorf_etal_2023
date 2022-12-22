######## Lipid Classes #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#        normalized, imputed, CV_filtered
#        Only Control (EV) and ash-2 RNAi(ASH) are plotted

# Approach: 1) Group lipids by classes (fatty acid information is discarded)
#           2) Fold change and stats 
#           3) Plots (TG, PE, PC, PI, Cer, DG). This code is repetitive. It only changes where marked.


# Goal: Plot lipid class changes of selected classes


# Output:  Pdf: 4_Class_TG.pdf
#          Csv: 4_Classes.csv   

#Packages
library(tidyverse)
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table_n<-df.impute_CVfiltered




######## 1) Class abundance all  ######################################################################################

class_table<-final_table_n%>%
  select(c(1,6:ncol(final_table_n))) #ignore fatty acid composition in columns

class_table<-aggregate(class_table[,-1], by = list(class_table$Class), FUN = sum) #aggregate all the classes


######## 2) Fold change and stats######################################################################################

#Mean, Coefficient of Variantion and Fold Change

class_table<-class_table %>% 
  mutate(ASH = apply(class_table %>% select(2:7), 1, mean ))%>% 
  mutate(EV = apply(class_table %>% select(8:13), 1, mean ))%>%
  mutate(PA = apply(class_table %>% select(14:19), 1, mean ))%>% 
  mutate(PRX = apply(class_table %>% select(20:25), 1, mean ))%>%
  mutate(cvASH = apply(class_table %>% select(2:7), 1, sd )/ASH)%>% 
  mutate(cvEV = apply(class_table %>% select(8:13), 1, sd )/EV)%>%
  mutate(cvPA = apply(class_table %>% select(14:19), 1, sd )/PA)%>% 
  mutate(cvPRX = apply(class_table %>% select(20:25), 1, sd )/PRX)%>%
  mutate(FC_ASHvsEV = (ASH/EV))%>%
  mutate(FC_PRXvsEV = (PRX/EV))%>%
  mutate(FC_PAvsEV = (PA/EV))%>%
  mutate(FC_PRXvsPA = (PRX/PA))


#Wilcoxon

for (i in 1:nrow(class_table)) {
  
  class_table$pval_ASHvsEV[i]<-wilcox.test(as.numeric(class_table[i,2:7]),as.numeric(class_table[i,8:13]))$p.value #comparing ASH vs EV
  class_table$pval_PAvsEV[i]<-wilcox.test(as.numeric(class_table[i,14:19]),as.numeric(class_table[i,8:13]))$p.value #comparing PRX-5_ASH vs EV
  class_table$pval_PRXvsEV[i]<-wilcox.test(as.numeric(class_table[i,20:25]),as.numeric(class_table[i,8:13]))$p.value #comparing PRX-5 vs EV
  class_table$pval_PRXvsPA[i]<-wilcox.test(as.numeric(class_table[i,20:25]),as.numeric(class_table[i,14:19]))$p.value #comparing PRX-5 vs PRX_ASH
}

#Adjust p-val

class_table<- class_table %>%
  mutate(FDR_ASHvsEV = p.adjust(pval_ASHvsEV,method = "BH"))%>%
  mutate(FDR_PAvsEV = p.adjust(pval_PAvsEV,method = "BH"))%>%
  mutate(FDR_PRXvsEV = p.adjust(pval_PRXvsEV,method = "BH"))%>%
  mutate(FDR_PRXvsPA = p.adjust(pval_PRXvsPA,method = "BH"))


#Add class names

colnames(class_table)[1]<-"Class"
rownames(class_table)<-class_table$Class

#Save
write.csv(class_table, '../Output_Data/4_Classes.csv')


######## Organize data for plot (only plotting EV, ash-2 RNAi)

tclass_table<-as.data.frame(t(class_table[,2:13])) %>%
  mutate(Group = NA)

tclass_table[c(1:6),ncol(tclass_table)]<- "ASH"

tclass_table[c(7:12),ncol(tclass_table)]<- "EV"

tclass_table$Group <- factor(tclass_table$Group , levels=c("EV", "ASH")) #order the groups

ptclass_table<-pivot_longer(tclass_table, cols = c(1:(ncol(tclass_table)-1)), names_to = "Var", values_to = "Val")



######## 3) Plots         ######################################################################################

######## Plot Triglycerides #################################################

plot<- filter(ptclass_table, grepl("TG",Var)) #changed class

Classes<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0, outlier.size=1.5)+
  geom_point(data=plot, aes(x=Group, y=Val), position=position_jitter(width=0.15)) +
  xlab("")+
  ylab("Triglyceride abundance (nmol/mg_prot)")+  #changed y-axis labelling
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

plot(Classes)

########Save Documents

ggsave("../Output_Figures/4_Class_TG.pdf", Classes, width=0.8, height=1.5, units="in", scale=3, useDingbats=FALSE) 



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
#   [1] data.table_1.13.6 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.5       purrr_0.3.4       readr_1.4.0       tidyr_1.1.3       tibble_3.1.1      ggplot2_3.3.3    
# [10] tidyverse_1.3.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8       cellranger_1.1.0 pillar_1.6.0     compiler_3.6.3   dbplyr_2.1.1     tools_3.6.3      digest_0.6.27    lubridate_1.7.10 jsonlite_1.7.2  
# [10] lifecycle_1.0.0  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.10     reprex_2.0.0     cli_2.5.0        DBI_1.1.0        rstudioapi_0.13  haven_2.3.1     
# [19] xml2_1.3.2       withr_2.4.2      httr_1.4.2       fs_1.5.0         generics_0.1.0   vctrs_0.3.7      hms_1.0.0        grid_3.6.3       tidyselect_1.1.0
# [28] glue_1.4.2       R6_2.5.0         fansi_0.4.1      readxl_1.3.1     farver_2.0.3     modelr_0.1.8     magrittr_2.0.1   MASS_7.3-53      backports_1.2.0 
# [37] scales_1.1.1     ellipsis_0.3.2   rvest_1.0.0      assertthat_0.2.1 colorspace_1.4-1 labeling_0.4.2   utf8_1.1.4       stringi_1.5.3    munsell_0.5.0   
# [46] broom_0.7.6      crayon_1.4.1
