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
  summarise(n=n(),mad=mad(tpm),med=median(tpm),q25=quantile(tpm,0.25),q75=quantile(tpm,0.75))

enh_peak_summary_tbl %>% 
  arrange(desc(med)) %>% 
  mutate(rank=dense_rank(med)) %>% 
  ggplot(.,aes(x=q25,xend=q75,y=rank,yend=rank))+
  geom_segment(size=0.05)+
  geom_point(aes(med,rank),color="red",size=0.05)+
  scale_x_log10()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/ENH_peak_IQR_segment.png")

enh_peak_summary_tbl %>% ggplot(.,aes(n))+geom_histogram()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/ENH_peak_nsample_hist.png")

enh_peak_summary_tbl %>% ggplot(.,aes(med,mad/med,color=n))+geom_point()+scale_x_log10()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/ENH_mad_vs_nsample_scatter.png")

enh_peak_summary_tbl<-cage_enh_long_tbl %>% 
  filter(tpm>0) %>% 
  group_by(Id) %>% 
  mutate(rank=min_rank(tpm)) %>% 
  summarise(gini=1 - (2/(n()-1))*(n()-(sum(rank*tpm))/(sum(tpm))),mad=mad(tpm),med=median(tpm),n=n())

enh_peak_summary_tbl %>% 
  filter(n>2) %>% 
  ggplot(.,aes(mad/med,gini))+
  geom_point(alpha=0.05)
