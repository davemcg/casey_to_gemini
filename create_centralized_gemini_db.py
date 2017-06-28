#!/usr/local/Anaconda/envs/py3.4/bin/python3

import argparse
import subprocess
import gzip
import glob
import time 

parser = argparse.ArgumentParser(description= 'Takes in list of vcf.gz (must be tabix\'ed) files from casey_xlsx_to_hgvs.py and sort_bgzip_tabix.sh and changes the sample name to something unique via a temporary vcf, then merges all of the files into a single vcf, and finally creates a gemini database from the master vcf.')

parser.add_argument('files',  help = 'vcf.gz files, comma separated. * allowed if double escaped')

args = parser.parse_args()

vcf_files = args.files
vcf_files = vcf_files.split(',')
new_vcf_file_names = []

# Opens each file
# creates a TEMP file in /scratch/mcgaugheyd 
# then writes the vcf to it with a unique sample name derived from the file name
# bgzips and tabix indexes each new TEMP file
for one_file_name in vcf_files:
	one_file = gzip.open(one_file_name, 'r')
	prefix = one_file_name.split('.vcf.gz')[0]
	temp_file = open('/scratch/mcgaugheyd/' + prefix + '.TEMP.vcf', 'w')
	vcf_data = one_file.read().decode('utf-8')
	vcf_data = vcf_data.replace('SAMPLE',prefix)
	temp_file.write(vcf_data)
	one_file.close()
	temp_file.close()

	subprocess.call('bgzip -f /scratch/mcgaugheyd/' + prefix + '.TEMP.vcf', shell = True)
	subprocess.call('tabix -f -p vcf /scratch/mcgaugheyd/' + prefix + '.TEMP.vcf.gz', shell = True)
	new_vcf_file_names.append('/scratch/mcgaugheyd/' + prefix + '.TEMP.vcf.gz')

# Now merges all the vcfs into one vcf
subprocess.call('module load vcftools', shell = True)
temp_master_vcf_name = '/scratch/mcgaugheyd/TEMP_casey_VCFs_' + time.time() + '.vcf'
subprocess.call('vcf-merge ' + ' '.join(new_vcf_file_names) + ' > ' + temp_master_vcf_name)

# Sort, bgzip, tabix
subprocess.call('/home/mcgaugheyd/git/casey_to_gemini/sort_bgzip_tabix.sh ' +t emp_master_vcf_name)
	






