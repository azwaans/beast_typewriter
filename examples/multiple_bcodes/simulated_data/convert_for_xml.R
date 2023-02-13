## ---------------------------
##
## Script name: nexus to txt
##
## Purpose of script: Convert the nexus format into manageable xml format
##
## Author: Antoine Zwaans
##
## Date Created: 2023-01-24
##
## Copyright (c) Antoine Zwaans, 2023
## Email: antoine.zwaans@bsse.ethz.ch
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac 

setwd("/Users/azwaans/all_beasts/beast_typewriter/examples/multiple_bcodes/simulated_data/")    # Antoine's working directory 

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
install.packages("ape")
library(ape)
library(stringr)


sink("alignment_gamma_1_20_bcodes.txt")

for(i in 1:10){
  con = file(paste0("small.clock0.0343711018512661.seed",i,".alignment.nexus"), "r")
  cat(paste0("<data  id=\"data_clock_0034_",i,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"14\"/>"))
  cat("\n")
  start_converting <- FALSE
  while ( TRUE ) {
    line = readLines(con, n = 1)
    
    if ( length(line) == 0 ) {
      break
    }
    if(start_converting && (line != "end;")) {
      convert_into_xml_line(line)  
    }
    if(line == "\tmatrix "){
      start_converting <- TRUE
    }
    
  }
  
  close(con)
  cat("\n")
  cat("</data>")
  cat("\n")
  cat("\n")
}

for(i in 1:10){
  
  con = file(paste0("small.clock0.1656288981487339.seed",i,".alignment.nexus"), "r")
  cat(paste0("<data  id=\"data_clock_0165_",i,"\" spec=\"Alignment\" name=\"alignment\" >
            <userDataType spec=\"beast.evolution.datatype.ScarData\" nrOfStates=\"14\"/>"))
  cat("\n")
  start_converting <- FALSE
  while ( TRUE ) {
    line = readLines(con, n = 1)
    
    if ( length(line) == 0 ) {
      break
    }
    if(start_converting && (line != "end;")) {
      convert_into_xml_line(line)  
    }
    if(line == "\tmatrix "){
      start_converting <- TRUE
    }
    
  }
  
  close(con)
  cat("\n")
  cat("</data>")
  cat("\n")
  cat("\n")
}


convert_into_xml_line <- function(line) {
  split_line <- str_split(line, " ")
  
  taxon_string <- gsub("\t","",split_line[[1]][1])
  sequence_string <- gsub(";", "",split_line[[1]][2])
  cat("
            ")
  cat((paste0("<sequence ","spec=\"Sequence\" taxon=\"",taxon_string,"\"  value=\"",sequence_string,"\"/>")))

}
