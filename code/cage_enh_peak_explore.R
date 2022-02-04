library(tidyverse)
library(vroom)
library(seriation)
library(viridisLite)
base::load("./data/sample_cell_line_ID_tbl.Rda")

cage_enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",comment = '#',col_select = contains(c("Id",sample_ID_tbl$sample.ID)))

cell_cor_mat<-cor(cage_enh_tbl[,-1],method = "spearman")

d<-as.dist(1-cell_cor_mat)
o<-seriate(d,method="HC")
image(cell_cor_mat[get_order(o),get_order(o)],col=viridis(100))
hist(as.numeric(cell_cor_mat))

cage_enh_long_tbl<-cage_enh_tbl %>% pivot_longer(cols=!Id, names_to = "sample.ID", values_to = "tpm")

enh_peak_summary_tbl<-cage_enh_long_tbl %>% 
  filter(tpm>0) %>% 
  group_by(Id) %>% 
  summarise(n=n(),mad=mad(tpm),med=median(tpm))

enh_peak_summary_tbl %>% ggplot(.,aes(n))+geom_density()

enh_peak_summary_tbl %>% ggplot(.,aes(med,mad/med,color=n))+geom_point()+scale_x_log10()


