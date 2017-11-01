#!/usr/local/apps/R/gcc_4.9.1/3.4.0/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(rmarkdown)
render(args[1],
		output_file = args[6],
		output_format = 'html_document',
        params = list(
        	master_db = args[2],
            patient_vcf = args[3],
			failed = args[4],
			patient_ID = args[5]),
		encoding = 'utf-8')
