library(renv)
renv::init()

renv::install("tidyverse")
renv::install("vroom")
renv::install("UpSetR")
renv::install("igraph")
renv::install("seriation")
renv::install("readxl")
renv::install("furrr")

renv::install("bioc::GenomicRanges")
renv::snapshot()
