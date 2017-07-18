#!/bin/bash

module load R/3.4.0_gcc-6.2.0
module load gemini/0.19.0

master_gemini_db=$1
sample_ID=$2


base_path='/data/mcgaugheyd/projects/nei/hufnagel/casey_panels/'
db=$base_path'gemini_db/'$master_gemini_db
patient_vcf=$base_path'vcf/'$sample_ID'.vt.vcf.gz'
failed=$base_path'FAILED/'$sample_ID'.FAILED.dat'
report_name=$base_path'reports/'$sample_ID'.report.html'

echo $report_name
Rscript ~/git/casey_to_gemini/src/write_report.R ~/git/casey_to_gemini/src/write_report.Rmd $db $patient_vcf $failed $sample_ID $report_name
