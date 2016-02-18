#!/usr/local/bin/python3

import pandas
import math
import os
import subprocess
import re

"""
Reads in xlsx file, extract hgvs, and builds ped
"""

# Reads all files in directory
all_files = os.listdir()
# only input patient xlsx files and builds a ped file (doesn't do gender)
xlsx = []
ped = []
ped.append('#Family SampleID PaternalID MaternalID Gender Phenotype Disease ClinicalID')
for i in all_files:
	if re.search('^\d+.*xlsx$',str(i)):
		xlsx.append(i)
		clinicalID = i.split()[0]	
		fname = i.split()[1]		
		lname = i.split()[2]
		name = fname + "_" + lname
		disease = i.split()[3]
		if "_" in disease:
			disease = disease.split("_")[0]
		ped_line = name + " " + name + " 0 " + "0 " + "0 " + "1 " + disease + " " + clinicalID
		ped.append(ped_line)



# loops through xlsx list, opens xlsx file to extract hgvs, builds dictionary of 
# hgvs as key and clinicalID_Zygosity as value

hgvs_dict = {}
for i in xlsx:
	# read first xlsx sheet into panda data structure (skipping the first 4 rows)
	data = pandas.read_excel(i, sheetname=0, skiprows=4)
	#just keep gene and hgvs coding info
	data = data[['Chromosome','HGVSGenomic','Zygosity']]
	# drop blanks (in any column)
	data = data.dropna()
	# drop anything that doesn't have a HGVS in it
	cond = ~data['HGVSGenomic'].str.contains('g.')
	data = data.drop(data[cond].index.values)
	# reprint with gene name
	data['HGVS'] = data["Chromosome"].map(str) + ":" + data['HGVSGenomic']
	
	
	# build dictionary
	clinicalID = i.split()[0]
	for index,row in data.iterrows():
		key = row['HGVS']
		if key in hgvs_dict:
			old_values = hgvs_dict[key]
			new_values = clinicalID + "_" + row['Zygosity']
			new_values = old_values + "," + new_values
			hgvs_dict[key] = new_values
			
		else:
			hgvs_dict[key] = clinicalID + "_" + row['Zygosity']
	
# print hgvs to file for vep
file=open("vep_hgvs_input.txt",'w')
for k,v in hgvs_dict.items():
	file.write(k)
	file.write('\n')

file.close()
	
	
######
# Run VEP. By default it left aligns the vcf
######
vep_query = 'perl /Applications/variant_effect_predictor/variant_effect_predictor.pl ' + \
		'-i vep_hgvs_input.txt ' + \
		'--vcf -o ' + 'vep_output.vcf ' + \
		'--assembly GRCh37 ' + \
		'--database' + ' --port 3337 ' + \
		'--plugin Grantham ' + \
		'--total_length ' + \
		'--hgvs ' + \
		'--sift b ' + \
    	'--polyphen b ' + \
    	'--symbol ' + \
		'--numbers ' + \
		'--biotype ' + \
		'--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Grantham,HGVSc,HGVSp ' + \
		'--force_overwrite' 

print("Running VEP online. Expect to wait ~10 minutes")
os.system(vep_query)
print("VEP query done!")

# import in vep output
vcf = open('vep_output.vcf','r')
vcf = vcf.readlines()
# rewrite vcf with sample genotypes
file = open('gemini.casey.vcf','w')
for line in vcf:
	if line[0:2] == '##':
		print(line[:-1])
	elif line[0:2] == '#C':
		
	

		
"""
Run VEP with hgvs, output vcf
"""

"""
Merge all vcfs into one vcf
"""

"""
Run VEP and annotate for gemini
"""

"""
Create gemini db, use ped
"""
		




