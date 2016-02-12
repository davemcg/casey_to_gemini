library(readxl)
library(dplyr)
setwd("~/Documents/PROJECTS/brooks/casey/")

casey_files <- list.files()

process_files <- function(file_name){
  data <- read_excel(file_name, skip=4)
  ID <- strsplit(file_name,"\s")[[1]][1]
  fname <- strsplit(file_name,"\s")[[1]][2]
  lname <- strsplit(file_name,"\s")[[1]][3]
  name <- paste(fname, lname)
  data$ID <- ID
  data$name <- name
  output <- data %>% select(ID,name,Gene,Chromosome,HGVSGenomic,HGVSCoding,Panel)
  return(data)
  
}


