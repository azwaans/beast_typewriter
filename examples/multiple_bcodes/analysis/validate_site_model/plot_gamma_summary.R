## ---------------------------
##
## Script name: Plot summaries from GammaCategories Site model runs
##
## Purpose of script:
##
## Author: Antoine Zwaans
##
## Date Created: 2022-12-21
##
## Copyright (c) Antoine Zwaans, 2022
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes: Based on S.Seidel scripts 
##
##
## ---------------------------

## set working directory for Mac

setwd("Users/azwaans/all_beasts/beast_typewriter/examples/multiple_bcodes/analysis/validate_site_model")
## ---------------------------


## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tracerer)
library(HDInterval)
library(ggplot2)

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R")

## ---------------------------


log_file = "typewriter_model_multiple_alignments_shared_insertRates_fromGamma1_small_truth.127.log"

log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)

data_gammaShape <- data.frame(log_data_wo_burnin$gammaShape)
colnames(data_gammaShape) <- "gammaShape"
data_gammaShape_all <- data.frame(name = as.factor(rep(1,length(data_gammaShape$gammaShape))),value = data_gammaShape)


for(i in 2:4) {
  
  log_file = paste0("typewriter_model_multiple_alignments_shared_insertRates_fromGamma",i,"_small_truth.127.log")
  log_data = parse_beast_tracelog_file(paste0(log_dir, log_file))
  log_data_wo_burnin = remove_burn_ins(log_data, burn_in_fraction = 0.1)
  
  data_gammaShape <- data.frame(log_data_wo_burnin$gammaShape)
  colnames(data_gammaShape) <- "gammaShape"
  data_gammaShape_all <- rbind(data_gammaShape_all,data.frame(name = as.factor(rep(i,length(data_gammaShape$gammaShape))),value = data_gammaShape))

}


g_shapes = ggplot(data_gammaShape_all, aes(x=name,y=gammaShape,fill=name)) +
  geom_violin() + ylab("GammaShape Estimates") + xlab("Truth") + ggtitle("GammaShape estimates - 4 categories")


true_shapes <- data.frame(name = c("1","2","3","4"),value=c(1.0,2.0,3.0,4.0))

g_shapes = g_shapes + geom_point(data=true_shapes,aes(x=name,y=value))

 
