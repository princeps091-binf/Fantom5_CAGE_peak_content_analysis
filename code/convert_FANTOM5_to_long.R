library(tidyverse)
library(vroom)

fantom5_tss_tbl<-"~/data_transfer/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt"

cage_tbl<-vroom(fantom5_tss_tbl,comment = '#')
col_ID<-str_split_fixed(colnames(cage_tbl)[-1],pattern = "\\.",4)
sample_ID<-col_ID[,3][-(1:6)]
colnames(cage_tbl)[-c(1:7)]<-sample_ID
cage_long_tbl<-cage_tbl %>% 
  dplyr::select(-c(2:7)) %>% 
  pivot_longer(!`00Annotation`,names_to="sample.ID",values_to="tpm")

save(cage_long_tbl,file="~/data_transfer/hg19.cage_peak_phase1and2combined_tpm_long.Rda")
