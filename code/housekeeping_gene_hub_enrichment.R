library(tidyverse)
library(vroom)
library(GenomicRanges)
library(readxl)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
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
fantom5_supp<-'~/Documents/multires_bhicect/data/epi_data/CAGE/41586_2014_BFnature13182_MOESM86_ESM.xlsx'
fantom5_tss_tbl<-"~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_long.Rda"
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"
spec_res_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
#------------------------------------------
H1<-c("CNhs14067","CNhs14068","CNhs13964")
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')
#------------------------------------------
cage_tbl<-input_tbl_fn(fantom5_tss_tbl)

housekeeping_peak_tbl<-read_excel(fantom5_supp,sheet = 5)
housekeeping_peak_set<-housekeeping_peak_tbl$`DPI peak`
cage_tbl %>% 
#  filter(`00Annotation` %in% housekeeping_peak_set) %>% 
  group_by(`00Annotation`) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n))+geom_density()

cell_cage_tbl<-cage_tbl %>% 
  filter(sample.ID %in% GM12878) %>% 
  group_by(`00Annotation`) %>% 
  summarise(m=mean(tpm)) %>% 
  filter(!(grepl("STAT:",`00Annotation`))) %>% 
  mutate(house=ifelse(`00Annotation` %in% housekeeping_peak_set,"house","other"))
cell_cage_tbl %>% 
  ggplot(.,aes(m,color=house))+
  geom_density()+
  scale_x_log10()

hub_tbl<-input_tbl_fn(hub_file) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])
hub_tbl<-do.call(bind_rows,map(unique(hub_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_folder,chromo,"_spec_res.Rda"))
  tmp_tbl<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,as.numeric)) 
  
})) %>% 
  dplyr::select(chr,node,res,bins)

plan(multisession,workers=4)
hub_tbl<-hub_tbl %>% 
  mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    inter_cl_Grange<-   GRanges(seqnames=chr,
                                ranges = IRanges(start=as.numeric(bins),
                                                 end=as.numeric(bins) + res_num[res]-1
                                ))
    inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
    return(inter_cl_Grange)
    
  }))
plan(sequential)

tmp_l<-hub_tbl %>% 
  #  filter(res == "100kb") %>% 
  #  filter(res%in% c("1Mb","500kb","100kb")) %>% 
  #  filter(res%in% c("10kb","50kb","5kb")) %>% 
  dplyr::select(GRange) %>% as.list

cl_GRange<-IRanges::reduce(do.call("c",tmp_l$GRange))

cell_cage_tbl<-cell_cage_tbl %>% 
  mutate(chr=str_split_fixed(`00Annotation`,":|\\.\\.|,",4)[,1],
         start=as.numeric(str_split_fixed(`00Annotation`,":|\\.\\.|,",4)[,2]),
         end=as.numeric(str_split_fixed(`00Annotation`,":|\\.\\.|,",4)[,3]))
cell_cage_Grange<-   GRanges(seqnames=cell_cage_tbl$chr,
                            ranges = IRanges(start=cell_cage_tbl$start,
                                             end=cell_cage_tbl$end
                            ))
mcols(cell_cage_Grange)<-tibble(m=cell_cage_tbl$m,house=cell_cage_tbl$house)

table(mcols(cell_cage_Grange)$house[unique(subjectHits(findOverlaps(cl_GRange,cell_cage_Grange)))])/length(unique(subjectHits(findOverlaps(cl_GRange,cell_cage_Grange))))
table(mcols(cell_cage_Grange)$house)/length(cell_cage_Grange)
obs<-as.numeric(table(mcols(cell_cage_Grange)$house[unique(subjectHits(findOverlaps(cl_GRange,cell_cage_Grange)))]))
exp_p<-as.numeric(table(mcols(cell_cage_Grange)$house)/length(cell_cage_Grange))
tmp_chi<-chisq.test(obs,p=exp_p)
tmp_chi$observed-tmp_chi$expected
