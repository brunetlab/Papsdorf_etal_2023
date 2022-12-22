######## SFA MUFA PUFA #################################################

# Data:2_Lipids_Imputed_CVfiltered (Output_Data)
#      Normalized, imputed, CV_filtered
#      Only the samples EV (EV) and ash-2 RNAi(ASH) are plotted

# Approach: 1) Fatty acid acyl chain abundance in all lipids (individual fatty acid concentration matches the corresponding complex lipid). Odd chain fatty acids and fatty acids with d, O and t modifications are not analyzed.
#           2) Categorize fatty acyl chain unsaturations into SFA, MUFA and PUFA
#           3) Calculate MUFA/PUFA
#           4) Fold change and stats 
#           5) Plot SFA, MUFA, PUFA, MUFAtoPUFA.  


# Goal: Plot SFA, MUFAs and PUFAs


# Output:  Pdf: 5_SFA_MUFA_PUFA_All_Lipids.pdf
#          Csv: 5_SFA_MUFA_PUFA_All_Lipids.csv   



#Packages
library(tidyverse)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))# Set wd to the current file location
rm(list=ls())

load(file = "../Output_Data/2_Lipids_Imputed_CVfiltered.Rdata")
final_table<-df.impute_CVfiltered


######## 1) Fatty acyl chain abundance in all lipids ###############################################################################

######## Fatty acid count. List with individual fatty acyl side chain listed with the corresponding lipid concentration in one table

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


#Table with the number of unsaturations (ether fatty acids are kept)

uns_FA_table<-uns_FA_table %>%
  mutate(FA = str_replace_all(FA, "e", "")) %>% #replace all the "e"s with spaces
  mutate(Uns = sub(".*:", "", FA)) %>% #delete the info before the ":"
  select(Uns,2:25)


#Aggregate

uns_FA_table<-aggregate(uns_FA_table[,-1], by = list(uns_FA_table$Uns), FUN = sum) #summarize intensity of the same fatty acids
uns_FA_table<-uns_FA_table[-1,]#delete the first row with the empty FA loadings 
uns_FA_table <-as.data.frame(lapply(uns_FA_table,as.numeric)) #make numeric



######## 2) Categorize unsaturations into SFA, MUFA and PUFA ##########################################################################

SMP<-t(uns_FA_table) #transpose
SMP<-as.data.frame(SMP)
colnames(SMP)<-SMP[1,] #set numer of unsaturations as colnames
SMP<-SMP %>% rename(SFA = 1, MUFA = 2) #rename SFA and MUFA
SMP<-SMP[-1,]
cols = 3:ncol(SMP) #set which columns are PUFAs
SMP[1:5, cols] # check
SMP$PUFA = rowSums(SMP[ , cols]) #summarize all >=2 to PUFA
head(SMP)
SMP = SMP[, c("SFA", "MUFA", "PUFA")]



######## 3) Calculate MUFA to PUFA ratio ##############################################################################################

SMP$MUFAtoPUFA = SMP[, c("MUFA")]/SMP[, c("PUFA")]
SMP_table<-t(SMP) #transpose
SMP_table<-as.data.frame(SMP_table)



######## 4) Fold change and stats #####################################################################################################

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

SMP_table$SMP = "SMP" 
SMP_table <- SMP_table %>%
  select(SMP, everything()) #move SMP to the front
SMP_table$SMP <-rownames(SMP_table)


#Save
write.csv(SMP_table, '../Output_Data/5_SFA_MUFA_PUFA_All_Lipids.csv')


########Organize data for plot. Plot only EV and ash-2 RNAi

tSMP_table<-as.data.frame(t(SMP_table[,2:13])) %>%
  mutate(Group = NA)

tSMP_table[c(1:6),ncol(tSMP_table)]<- "ASH"

tSMP_table[c(7:12),ncol(tSMP_table)]<- "EV"

tSMP_table$Group <- factor(tSMP_table$Group , levels=c("EV", "ASH")) #order the groups

ptSMP_table<-pivot_longer(tSMP_table, cols = c(1:(ncol(tSMP_table)-1)), names_to = "Var", values_to = "Val")

ptSMP_table$Var <- factor(ptSMP_table$Var , levels=c ("SFA", "MUFA", "PUFA", "MUFAtoPUFA")) #order values



######## 5) Plots         ###########################################################################################################


######## Plot SFA         #################################################

plot<- filter(ptSMP_table, grepl("SFA",Var)) #subset SFA


######## plot SFA

SFA<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0, outlier.size=1.5)+
  geom_point(data=plot, aes(x=Group, y=Val), position=position_jitter(width=0.15)) +
  xlab("")+
  ylab("SFA abundance all lipids (nmol/mg_prot)")+  #changed 
  ylim(70,200)+ #to have always the same scale across SFA MUFA PUFA
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

plot(SFA)


######## Plot MUFA         #################################################

plot<- filter(ptSMP_table, grepl("MUFA",Var)) #subset MUFA
plot<- filter(plot, !grepl("MUFAtoPUFA",Var)) 


######## plot MUFA 

MUFA<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0, outlier.size=1.5)+
  geom_point(data=plot, aes(x=Group, y=Val), position=position_jitter(width=0.15)) +
  xlab("")+
  ylab("MUFA abundance all lipids(nmol/mg_prot)")+  #changed
  ylim(70,200)+ #to have always the same axis across SFA MUFA PUFA
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

plot(MUFA)



######## Plot PUFA         #################################################

plot<- filter(ptSMP_table, grepl("PUFA",Var)) #subset PUFA
plot<- filter(plot, !grepl("MUFAtoPUFA",Var)) 


######## plot PUFA

PUFA<-ggplot(plot,aes(x = Group, y = Val, fill = Group)) +
  geom_boxplot(outlier.alpha = 0, outlier.size=1.5)+
  geom_point(data=plot, aes(x=Group, y=Val), position=position_jitter(width=0.15)) +
  xlab("") +
  ylab("PUFA abundance all lipids (ng/mg_prot)")+  #changed
  ylim(70,200)+ #to have always the same axis across SFA MUFA PUFA
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
                               "ASH" = "red1",
                               "PRX" = "blue"))

plot(PUFA)


######## Plot all together  #################################################

SFA_MUFA_PUFA<-ggarrange(SFA, MUFA, PUFA, 
                    ncol = 3, nrow = 1)

plot(SFA_MUFA_PUFA)

ggsave("../Output_Figures/5_SFA-MUFA-PUFA_All_Lipids.pdf", SFA_MUFA_PUFA, width=3.2, height=2, units="in", scale=3)





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
# > 
