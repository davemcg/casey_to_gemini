#!/usr/local/apps/R/gcc_6.2.0/3.4.0/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

rmarkdown::render(args[1],
				  output_file = args[6],
                  params = list(
                      master_db = args[2],
                      patient_vcf = args[3],
					  failed = args[4],
					  patient_ID = args[5])
                  )
