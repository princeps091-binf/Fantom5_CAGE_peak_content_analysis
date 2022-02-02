library(tidyverse)
library(vroom)
library(igraph)
library(seriation)

cage_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#',col_select = contains(c("Annotation","cell%20line")))
cell_line<-colnames(cage_tbl)
sample_id<-unlist(lapply(strsplit(cell_line,split='\\.'),'[',3))
long_transform_tbl<-do.call(bind_rows,lapply(2:ncol(cage_tbl),function(i){
  
  cage_tbl %>% dplyr::select(1,i) %>%
    slice(-c(1,2)) %>% 
    mutate(value=.data[[cell_line[i]]]) %>%
    dplyr::select(-2) %>% 
    filter(value > 0) %>% 
    mutate(sample.ID=sample_id[i])
  
}))


sample_summary_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(sample.ID) %>% 
  summarise(n=n(),IQR=IQR(value),mad=mad(value),med=median(value),CV=sd(value)/mean(value),CV2=mad(value)/median(value))


sample_summary_tbl %>% 
  ggplot(.,aes(med,IQR))+
  geom_point()+scale_x_log10()

sample_summary_tbl %>% 
  ggplot(.,aes(n,IQR))+
  geom_point()+scale_x_log10()

sample_summary_tbl %>% 
  ggplot(.,aes(n,CV2))+
  geom_point(alpha=0.3)

edge_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(sample.ID) %>% 
  mutate(weight=1/dense_rank(desc(value))) %>% 
  dplyr::select(peak.ID,sample.ID,weight)
ego_tbl<-edge_tbl %>% ungroup() %>% 
  distinct(peak.ID) %>% mutate(ego.idx=1:n())
alter_tbl<-edge_tbl %>% ungroup %>% 
  distinct(sample.ID) %>% mutate(alter.idx=1:n())
edge_tbl<-edge_tbl %>% ungroup %>% 
  left_join(.,ego_tbl) %>% 
  left_join(.,alter_tbl)

FANTOM_mat<-Matrix::sparseMatrix(i=edge_tbl$ego.idx,j=edge_tbl$alter.idx,x=log10(edge_tbl$weight))
png("./img/FANTOM_sparsity.png")
image(as.matrix(FANTOM_mat))
dev.off()

o <- c(

    seriate(dist(FANTOM_mat, "minkowski", p = 1), method = "TSP"),

    seriate(dist(Matrix::t(FANTOM_mat), "minkowski", p = 1), method = "TSP")
)
