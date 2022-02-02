library(renv)
renv::init()

renv::install("tidyverse")
renv::install("vroom")

renv::snapshot()
