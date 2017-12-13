#!/usr/local/apps/R/gcc_4.9.1/3.4.0/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(rmarkdown)
render(input = args[1],
		output_file = args[5],
		output_format = 'html_document',
        params = list(
          gemini_return = args[2],
          patient_MVL_excel = args[3],
			patient_ID = args[4]),
		encoding = 'utf-8')
