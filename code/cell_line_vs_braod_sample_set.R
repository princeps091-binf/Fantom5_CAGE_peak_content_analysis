library(tidyverse)
library(vroom)
#------------------------------------------
## Utils. Fn
input_tbl_fn<-function(file){
in_tbl<-  get(load(file))
tmp_obj<-names(mget(load(file)))
rm(list=tmp_obj)
rm(tmp_obj)
return(in_tbl)
}
#------------------------------------------
small_set_file<-'./data/long_form_cage_tbl.Rda'
big_set_file<-'./data/long_form_big_set_cage_tbl.Rda'
small_set_tbl<-input_tbl_fn(small_set_file) %>% group_by(`00Annotation`) %>% summarise(n=n())
big_set_tbl<-input_tbl_fn(big_set_file) %>% group_by(peak.ID) %>% summarise(n=n())
gg_scat<-big_set_tbl %>%
full_join(.,small_set_tbl %>% dplyr::rename(peak.ID=`00Annotation`)
,by=c("peak.ID")) %>%
ggplot(.,aes(n.x,n.y))+geom_point(alpha=0.05)+xlab("active samples in broad sample set")+ylab("active sample in cell-line sample set")
ggsave("~/Documents/multires_bhicect/weeklies/weekly51/img/broad_vs_line_scatter.png")
