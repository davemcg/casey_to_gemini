#!/usr/local/Anaconda/envs/py3.4/bin/python3

import argparse
import subprocess
import gzip
import glob
import time 

parser = argparse.ArgumentParser(description= 'Takes in list of vcf.gz (must be tabix\'ed) files from casey_xlsx_to_hgvs.py and sort_bgzip_tabix.sh and changes the sample name to something unique via a temporary vcf, then merges all of the files into a single vcf, and finally creates a gemini database from the master vcf.')

parser.add_argument('files',  help = 'vcf.gz files, comma separated. * allowed if double escaped')
parser.add_argument('gemini_db_name', help = 'name for Gemini database')
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
time_stamp = str(time.time())
temp_master_vcf_name = '/scratch/mcgaugheyd/TEMP_casey_VCFs_' + time_stamp + '.vcf'
subprocess.call('vcf-merge ' + ' '.join(new_vcf_file_names) + ' > ' + temp_master_vcf_name, shell = True)

# Sort, bgzip, tabix
subprocess.call('/home/mcgaugheyd/git/casey_to_gemini/sort_bgzip_tabix.sh ' + temp_master_vcf_name, shell = True)

# create fake ped 
temp_ped_file = '/scratch/mcgaugheyd/TEMP_casey_' + time_stamp + '.ped'
temp_ped = open(temp_ped_file, 'w')
for sample in vcf_files:
	temp_ped.write(sample.split('.vcf.gz')[0] + ' ' + sample.split('.vcf.gz')[0] + ' 0 0 0 2\n')
temp_ped.close()

# move temp files to current dir (GATK_vcf_to_geminiDB.sh requires this)
subprocess.check_call(['cp', temp_ped_file, '.'])
subprocess.check_call(['cp', temp_master_vcf_name + '.gz', '.'])

# create db
subprocess.call('sbatch --partition=quick --cpus-per-task 8 /home/mcgaugheyd/git/variant_prioritization/GATK_vcf_to_geminiDB.sh ' + temp_master_vcf_name.split('/')[-1] + '.gz ' + temp_ped_file.split('/')[-1] + ' ' + args.gemini_db_name, shell = True)

# remove temp files
subprocess.call(['rm',temp_ped_file])
subprocess.call(['rm',temp_master_vcf_name + '.gz'])







