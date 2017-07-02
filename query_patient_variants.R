#!/usr/bin/env Rscript

# Takes in master gemini DB and patient-specific vcf (the VEP/vcfanno processed version)
# and returns a data frame with the patient variants

args = commandArgs(trailingOnly=TRUE)

gemini_db <- args[1]
file <- args[2]


library(data.table)
library(tidyverse)
#gemini_db <- 'test_40.db'
#file <- '17-1012_RD_v11_assay.VEP.GRCh37.anno.vcf.gz'
setwd('/data/mcgaugheyd/projects/nei/hufnagel/casey_panels/test2')

con <- file(file, open = 'r')
data <- readLines(con)
row_skip <- sum(sapply(data, function(x) ifelse(substr(x,1,2)=='##', 1,0)))

vcf <- fread(paste('zcat', file), skip=row_skip)
id_status <- vcf %>% select(`#CHROM`, POS, ID, REF, ALT, SAMPLE) %>% separate(SAMPLE, c('GT', 'CT', 'JC'), ':')

primary_IDs <- id_status %>% filter(JC=='Primary') %>% .[['ID']]
secondary_IDs <- id_status %>% filter(JC=='Secondary') %>% .[['ID']]

gemini_query <- function(gemini_db, vcf_id, ref, alt){
  # Queries gemini in two ways:
  ## With the HGVS as the search term
  ## If ref, alt given, then uses HGVS as a fuzzy match and ref alt to confirm
  base_query <- "gemini query --header -q 'select chrom, start, end, ref, alt, vcf_id, hgvsc, hgvsp, type, gene, impact, rs_ids, pfam_domain, clinvar_diseases, clinvar_ID, clinvar_pathogenic, clinvar_sig, hgmd_overlap, af_exac_afr, af_exac_all af_exac_amr, af_exac_eas, af_exac_nfe af_exac_oth, af_exac_sas, an_exac_all, exac_num_het, exac_num_hom_alt, mis_z, pli, gerp_elements, polyphen_pred, polyphen_score, sift_pred, sift_score, cadd_phred FROM variants "

  if (missing(ref) & missing(alt)){
       out <- system(paste0(base_query,
                "WHERE vcf_id IN (", 
                         paste0('\"',vcf_id, '\"', collapse=','), ") ' ", 
                         gemini_db), 
                  intern=TRUE) 
  }
  else {
    out <- system(paste0(base_query,
                         "WHERE vcf_id LIKE ", 
                         paste0('\"%',vcf_id, '%\"', collapse=','), 
                         " AND ref IS \"", ref, 
                         "\" AND alt IS \"", alt,"\"' ", 
                         gemini_db), 
                  intern=TRUE) 
  }
  return(out)
}

DF_maker <- function(df_vector, priority){
  # Converts character vector of gemini captured call into a data frame
  DF <- data.frame(do.call(rbind, strsplit(df_vector, '\t')))
  colnames(DF) <- as.character(unlist(DF[1,]))
  DF <- DF[-1,] 
  DF$Priority <- priority
  return(DF)
}
ID_to_DF <- function(IDs, priority){
  # Converts ID (HGVS) vector to data frame with status 
  id_DF <- data.frame(IDs)
  id_DF$Priority <- priority
  colnames(id_DF)[1] <- 'HGVS'
  return(id_DF)
}


primary_DF <- DF_maker(gemini_query(gemini_db, primary_IDs), 'Primary')
secondary_DF <- DF_maker(gemini_query(gemini_db, secondary_IDs), 'Secondary')
all_DF <- rbind(primary_DF, secondary_DF)
primary_IDs <- ID_to_DF(primary_IDs, 'Primary')
secondary_IDs <- ID_to_DF(secondary_IDs, 'Secondary')
all_IDs <- rbind(primary_IDs, secondary_IDs)

#add failures to match ID/vcf_id 
HGVS = all_IDs$HGVS
if (length(HGVS[!HGVS %in% all_DF$vcf_id]) > 0){
  for (i in HGVS[!HGVS %in% all_DF$vcf_id]){
    ref <- id_status %>% filter(ID==i) %>% .[['REF']]
    alt <- id_status %>% filter(ID==i) %>% .[['ALT']]
    query_result <- gemini_query(gemini_db, i, ref, alt)
    if (length(query_result) != 2) {
      print(paste("ERROR: missing ", i)) 
      quit(status = 1)
    }
    priority <- all_IDs %>% filter(HGVS==i) %>% pull(Status)
    all_DF <- rbind(all_DF, DF_maker(gemini_query(gemini_db, i, ref, alt), priority) )
  }
}

all_DF

