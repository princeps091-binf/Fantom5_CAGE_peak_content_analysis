library(tidyverse)
library(vroom)
cage_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#',col_select = contains(c("Annotation","cell%20line")))
cell_line<-colnames(cage_tbl)
strsplit(unlist(lapply(strsplit(cell_line,split='\\.'),'[',2)),split="%20cell%20line%3a")
sample_id<-unlist(lapply(strsplit(cell_line,split='\\.'),'[',3))
prank_transform_tbl<-do.call(bind_rows,lapply(2:ncol(cage_tbl),function(i){
  
  cage_tbl %>% dplyr::select(1,i) %>%
    slice(-c(1,2)) %>% 
    filter(.data[[cell_line[i]]]>0) %>% 
    mutate(p=percent_rank(.data[[cell_line[i]]])) %>% 
    dplyr::select(-2) %>% 
    rename(value=p) %>% mutate(sample.ID=sample_id[i])
  
}))
