library(tidyverse)
library(vroom)
library(igraph)
library(seriation)
library(viridisLite)

cage_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#',col_select = contains(c("Annotation","cell%20line")))
cell_line<-colnames(cage_tbl)
sample_id<-unlist(lapply(strsplit(cell_line,split='\\.'),'[',3))
colnames(cage_tbl)<-c("peak.ID",sample_id[-1])
cage_tbl<-cage_tbl[-c(1,2),]
#Produce correlation matrix across cell-lines
cell_cor_mat<-cor(cage_tbl[,-1],method = "spearman")
save(cell_cor_mat,file = './data/CAGE_cormat.Rda')

base::load('./data/CAGE_cormat.Rda')
d<-as.dist(1-cell_cor_mat)
o<-seriate(d,method="HC")
image(cell_cor_mat[get_order(o),get_order(o)],col=viridis(100))

cell_cor_mat[cell_cor_mat<median(as.numeric(cell_cor_mat))]<-0

g_cell<-graph_from_adjacency_matrix(cell_cor_mat,weighted = T,diag = F,mode = "undirected")
#two samples isolated once we set correlation threshold
main_g_v<-names(which(components(g_cell)$membership==1))

main_sub_g<-induced_subgraph(g_cell,main_g_v)
louvain_sample_cluster<-cluster_louvain(main_sub_g)
table(cluster_louvain(main_sub_g)$membership)
grep(paste(V(main_sub_g)$name[which(cluster_louvain(main_sub_g)$membership==3)],collapse = "|"),cell_line,value=T)
##Build the corresponding clustering matrix
sample_comm<-unique(louvain_sample_cluster$membership)
comm_edge_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-V(main_sub_g)$name[which(cluster_louvain(main_sub_g)$membership==i)]
  expand_grid(ego=tmp_v,alter=tmp_v) %>% mutate(x=i)
}))
tmp_mat<-matrix(0,nrow = nrow(cell_cor_mat),ncol=ncol(cell_cor_mat),dimnames = dimnames(cell_cor_mat))
tmp_mat[as.matrix(comm_edge_tbl[,1:2])]<-comm_edge_tbl$x
image(tmp_mat[get_order(o),get_order(o)],col=plasma(4))

sample_cl_tbl<-do.call(bind_rows,lapply(sample_comm,function(i){
  tmp_v<-V(main_sub_g)$name[which(cluster_louvain(main_sub_g)$membership==i)]
  return(tibble(sample.ID=tmp_v,cl=i))
}))

sample_cl_tbl<-tibble(sample.ID=V(g_cell)$name) %>% 
  left_join(.,sample_cl_tbl)
#Pivot to more efficient long-format
long_transform_tbl<-do.call(bind_rows,lapply(2:ncol(cage_tbl),function(i){
  
  cage_tbl %>% dplyr::select(1,i) %>%
    slice(-c(1,2)) %>% 
    mutate(value=.data[[cell_line[i]]]) %>%
    dplyr::select(-2) %>% 
    filter(value > 0) %>% 
    mutate(sample.ID=sample_id[i])
  
}))

save(long_transform_tbl,file = './data/long_form_cage_tbl.Rda')
base::load('./data/long_form_cage_tbl.Rda')

long_transform_tbl<-long_transform_tbl %>% 
  left_join(.,sample_cl_tbl)

sample_ID_tbl<-long_transform_tbl %>% distinct(sample.ID)
save(sample_ID_tbl,file = './data/sample_cell_line_ID_tbl.Rda')


long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(!(is.na(cl))) %>% 
  ggplot(.,aes(value,color=as.factor(cl)))+geom_density()+scale_x_log10()

long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  ggplot(.,aes(value))+geom_density()+scale_x_log10()


long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(!(is.na(cl))) %>% 
  filter(peak.ID=="chr1:10003486..10003551,+") %>% 
  ggplot(.,aes(value,color=as.factor(cl)))+geom_density()


gg_peak_cl_med<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(!(is.na(cl))) %>% 
  group_by(peak.ID,cl) %>% 
  summarise(med=median(value),n=n()) %>%
  filter(n>10) %>% 
  ggplot(.,aes(med,peak.ID,color=as.factor(cl)))+geom_point()+scale_x_log10()+
  theme(axis.ticks.y=element_blank())

ggsave("./img/peak_cl_med.png",gg_peak_cl_med)

long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(!(is.na(cl))) %>% 
  group_by(peak.ID,cl) %>% 
  summarise(med=median(value),n=n()) %>%
  ungroup() %>% 
  group_by(peak.ID) %>% 
  summarise(sd=abs(diff(range(med))),m=mean(med),n=n()) %>%
  filter(n>1) %>% 
  mutate(cv=sd/m) %>% 
  ggplot(.,aes(cv,m))+geom_point()+scale_x_log10()+scale_y_log10()

sample_summary_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(sample.ID) %>% 
  summarise(n=n(),IQR=IQR(value),mad=mad(value),med=median(value),CV=sd(value)/mean(value),CV2=mad(value)/median(value))


peak_summary_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>%
#  filter(!(is.na(cl))) %>% 
  group_by(peak.ID) %>% 
  summarise(n=n(),mad=mad(value),med=median(value),q25=quantile(value,0.25),q75=quantile(value,0.75))

peak_summary_tbl %>% 
  arrange(desc(med)) %>% 
  mutate(rank=dense_rank(med)) %>% 
  ggplot(.,aes(x=q25,xend=q75,y=rank,yend=rank))+
  geom_segment(size=0.05)+
  geom_point(aes(med,rank),color="red",size=0.05)+
  scale_x_log10()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/TSS_peak_IQR_segment.png")

peak_summary_tbl %>% 
  ggplot(.,aes(n))+geom_histogram()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/TSS_peak_nsample_hist.png")

peak_summary_tbl %>% 
  ggplot(.,aes(mad/med))+geom_histogram()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/TSS_peak_npcv_hist.png")

peak_summary_tbl %>% ggplot(.,aes(med,mad/med,color=n))+geom_point()+scale_x_log10()
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/TSS_mad_vs_nsample_scatter.png")

peak_summary_tbl %>% mutate(nc=ifelse(n>300,"ubiquitous",ifelse(n<50,"specialised","intermediate"))) %>% 
  mutate(nc=fct_relevel(nc,c("specialised","intermediate","ubiquitous"))) %>% 
  ggplot(.,aes(med,mad/med,color=n))+geom_point(alpha=0.1)+scale_x_log10()+facet_grid(nc~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/TSS_mad_vs_nsample_ncated_scatter.png")

#------------------------------
## Lorenz curve
tmp_sample<-"CNhs12331"
gg_lorenz<-long_transform_tbl %>% 
  group_by(sample.ID) %>% 
  filter(n()>5) %>% 
  mutate(p=percent_rank(value)) %>% 
  arrange(p) %>% 
  mutate(lo=cumsum(value)) %>% 
  mutate(lo=lo/sum(value)) %>% 
  ggplot(.,aes(p,lo,group=sample.ID))+geom_line(size=0.1)+
  geom_abline(slope = 1,intercept = 0, color="red")
ggsave("./img/TSS_lorenz.png",gg_lorenz)  

tmp_peak<-"chr10:100174900..100174956,-"

gg_lorenz<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(peak.ID==tmp_peak) %>% 
  filter(n()>300) %>% 
  mutate(p=percent_rank(value)) %>% 
  arrange(p) %>% 
  mutate(lo=cumsum(value)) %>% 
  mutate(lo=lo/sum(value)) %>% 
  ggplot(.,aes(p,lo,group=peak.ID))+geom_line(size=0.1,alpha=0.1)+
  geom_abline(slope = 1,intercept = 0, color="red")
ggsave("./img/TSS_ubiquitous_lorenz.png",gg_lorenz)  
## compute the gini coef

peak_gini_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(peak.ID) %>% 
  mutate(rank=min_rank(value)) %>% 
  summarise(gini=1 - (2/(n()-1))*(n()-(sum(rank*value))/(sum(value))),n=n(),mad=mad(value),med=median(value)) 

peak_gini_tbl %>% 
  filter(is.na(gini)) %>% arrange(desc(n))

peak_gini_tbl %>% 
  ggplot(.,aes(gini))+geom_density()

peak_gini_tbl %>% 
  ggplot(.,aes(n,gini))+geom_point(alpha=0.01)

peak_gini_tbl %>% 
  filter(n>2) %>% 
  ggplot(.,aes(med,gini,col=n))+
  scale_x_log10()+
  geom_point(alpha=0.05)

peak_gini_tbl %>% 
  filter(n>2) %>% 
  ggplot(.,aes(mad/med,gini))+geom_point(alpha=0.05)

peak_gini_tbl%>% 
  filter(n>2) %>% 
  mutate(nc=ifelse(n>300,"ubiquitous",ifelse(n<50,"specialised","intermediate"))) %>% 
  mutate(nc=fct_relevel(nc,c("specialised","intermediate","ubiquitous"))) %>% 
  ggplot(.,aes(med,gini,col=n))+
  scale_x_log10()+
  geom_point(alpha=0.05)+
  facet_grid(nc~.)


peak_gini_tbl%>% 
  filter(n>2) %>% 
  mutate(nc=ifelse(n>300,"ubiquitous",ifelse(n<50,"specialised","intermediate"))) %>% 
  mutate(nc=fct_relevel(nc,c("specialised","intermediate","ubiquitous"))) %>% 
  ggplot(.,aes(mad/med,gini,col=n))+
  geom_point(alpha=0.05)+
  facet_grid(nc~.)

tmp_peak<-"chr1:100110807..100110818,+"
long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(peak.ID == tmp_peak) %>% 
  ggplot(.,aes(value))+geom_histogram()


sample_gini_tbl<-long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  group_by(sample.ID) %>% 
  mutate(rank=min_rank(value)) %>% 
  summarise(gini=1 - (2/(n()-1))*(n()-(sum(rank*value))/(sum(value))),n=n(),mad=mad(value),med=median(value)) 

sample_gini_tbl %>% 
  filter(is.na(gini)) %>% arrange(desc(n))

sample_gini_tbl %>% 
  ggplot(.,aes(gini))+geom_density()

sample_gini_tbl %>% 
  ggplot(.,aes(n,gini))+geom_point()
