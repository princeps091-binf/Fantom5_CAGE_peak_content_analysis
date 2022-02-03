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

long_transform_tbl %>% 
  dplyr::rename(peak.ID=`00Annotation`) %>% 
  filter(!(is.na(cl))) %>% 
  ggplot(.,aes(value,color=as.factor(cl)))+geom_density()+scale_x_log10()


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
  filter(!(is.na(cl))) %>% 
  group_by(peak.ID) %>% 
  summarise(n=n(),mad=mad(value),med=median(value),q25=quantile(value,0.25),q75=quantile(value,0.75))

peak_summary_tbl %>% ggplot(.,aes(med,mad/med,color=n))+geom_point()+scale_x_log10()

peak_summary_tbl %>% mutate(nc=ifelse(n>300,"ubiquitous",ifelse(n<10,"specialised","else"))) %>% 
  ggplot(.,aes(med,mad/med,color=n))+geom_point(alpha=0.1)+scale_x_log10()+facet_grid(nc~.)

peak_summary_tbl %>% 
  ggplot(.,aes(n))+geom_density()

peak_summary_tbl %>% 
  arrange(desc(med)) %>% 
  mutate(rank=dense_rank(med)) %>% 
  ggplot(.,aes(x=q25,xend=q75,y=rank,yend=rank))+
  geom_segment(size=0.05)+
  geom_point(aes(med,rank),color="red",size=0.1)+
scale_x_log10()

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
