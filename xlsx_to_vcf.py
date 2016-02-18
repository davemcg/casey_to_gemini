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



# loops through xlsx list, opens xlsx file to extract hgvs
for i in xlsx:
	# read first xlsx sheet into panda data structure (skipping the first 4 rows)
	data = pandas.read_excel(i, sheetname=0, skiprows=4)
	#just keep gene and hgvs coding info
	data = data[['Chromosome','HGVSGenomic']]
	# drop blanks
	data = data.dropna()
	# drop anything that doesn't have a HGVS in it
	cond = ~data['HGVSGenomic'].str.contains('g.')
	data = data.drop(data[cond].index.values)
	# reprint with gene name
	data['HGVS'] = data["Chromosome"].map(str) + ":" + data['HGVSGenomic']
	
	#write to temp file
	file = open('vep_temp.txt','w')
	for i in data['HGVS']:
		file.write(i)
		file.write('\n')
	file.close()
	######
	# Run VEP
	######
	clinicalID = i.split()[0]
	vcf_name = clinicalID + '.vcf '
	vep_query = 'perl /Applications/variant_effect_predictor/perl variant_effect_predictor.pl ' + \
		'-i vep_temp.txt ' + \
		'--vcf -o ' + vcf_name + \
		'--species human --assembly GRCh37 ' + \
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
		




