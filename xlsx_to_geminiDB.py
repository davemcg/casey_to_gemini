#!/usr/local/bin/python3

import pandas
import math
import os
import subprocess
import re
from datetime import datetime
from collections import defaultdict

##########################################################################################
#Reads in xlsx file, extract hgvs, and writes ped file
##########################################################################################

# Reads all files in directory
all_files = os.listdir()
# only input patient xlsx files and builds a ped file (doesn't do gender)
xlsx = []
ped = []
ped.append('#Family SampleID PaternalID MaternalID Gender Phenotype Disease ClinicalID')

#build list of all samples for use through script
all_names = []
for i in all_files:
	# ID file names that start with a digit and end in xlsx
	if re.search('^\d+.*xlsx$',str(i)):
		xlsx.append(i)
		clinicalID = i.split()[0]
		# extract first and last name from file name
		fname = i.split()[1]		
		lname = i.split()[2]
		name = fname + '_' + lname
		all_names.append(name)
		# and pull disease from file name
		disease = i.split()[3]
		# clean up a formatting thing where some diseases are appended with "_"
		if '_' in disease:
			disease = disease.split('_')[0]
		ped_line = name + ' ' + name + ' 0 ' + '0 ' + '0 ' + '1 ' + disease + ' ' + clinicalID
		ped.append(ped_line)

print('\n\n' + ','.join(all_names) + ' sample HGVS extracted and ped built\n')

file = open('casey.ped', 'w')
for line in ped:
	file.write(line)
	file.write('\n')

file.close()
print('casey.ped written\n')

##########################################################################################
# Builds sample dictionary with hgvs, clinicalID, and allele status
##########################################################################################

# loops through xlsx list, opens xlsx file to extract hgvs, builds dictionary of 
# hgvs as key and clinicalID_Zygosity as value
panels = []
hgvs_dict = defaultdict(dict)
for i in xlsx:
	# check to find how many rows to skip
	check = pandas.read_excel(i,sheetname=0)
	# pull first column and convert to list
	first_column = list(check[check.columns[0]])
	index = first_column.index('#ID') + 1
	# read first xlsx sheet into panda data structure (skipping the first 4 rows)
	data = pandas.read_excel(i, sheetname=0, skiprows=index)
	#just keep gene and hgvs coding info
	data = data[['HGVSCoding','Zygosity','TimesObservedPerPanel','Panel']]
	# drop blanks (in any column)
	data = data.dropna()
	# drop anything that doesn't have a HGVS in it
	cond = ~data['HGVSCoding'].str.contains('c.')
	data = data.drop(data[cond].index.values)
	# build defaultdictionary
	fname = i.split()[1]		
	lname = i.split()[2]
	name = fname + "_" + lname
	for index,row in data.iterrows():
		key = row['HGVSCoding']
		key1, key2 = key.split(':')
		key1 = key1.split('.')[0]
		if key1 == ' ' or key1 == '':
			print(i, key1, key2, "Missing gene!!!!!!\n\nNot being entered into db!!!!!")
			continue
		key = key1 + ':' + key2
		panels.append(row['Panel'])
		hgvs_dict[key][name] = row['Zygosity'] + '_' +  str(row['TimesObservedPerPanel']) + '_' + row['Panel']


print("\n\nHere are the panels used across the files: ", set(panels), "\n\n")

# print hgvs to file for counsyl hvgs (hgvsC_to_vcf-like.py)
file=open("hgvs_input.txt",'w')
for k,v in hgvs_dict.items():
	file.write(k)
	file.write('\n')

file.close()
##########################################################################################
# Convert to VEP friendly format and re-do hgvs_dict with chr_pos_ref_alt as key
##########################################################################################
os.system('~/git/casey_to_gemini/./hgvsC_to_vcf-like.py hgvs_input.txt')


##########################################################################################
# Run online VEP. By default it left aligns the vcf(?)
##########################################################################################

# will output vcf from the hgvs
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


##########################################################################################
# add conversion info from hgvsC_to_vcf-like.py into 



##########################################################################################
# rewrite vcf with sample genotypes
##########################################################################################

# import in vep output
vcf = open('vep_output.vcf','r')
vcf = vcf.readlines()

# new header lines
new_header = \
'##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n\
##contig=<ID=chrM,length=16571,assembly=hg19>\n\
##contig=<ID=chr1,length=249250621,assembly=hg19>\n\
##contig=<ID=chr2,length=243199373,assembly=hg19>\n\
##contig=<ID=chr3,length=198022430,assembly=hg19>\n\
##contig=<ID=chr4,length=191154276,assembly=hg19>\n\
##contig=<ID=chr5,length=180915260,assembly=hg19>\n\
##contig=<ID=chr6,length=171115067,assembly=hg19>\n\
##contig=<ID=chr7,length=159138663,assembly=hg19>\n\
##contig=<ID=chr8,length=146364022,assembly=hg19>\n\
##contig=<ID=chr9,length=141213431,assembly=hg19>\n\
##contig=<ID=chr10,length=135534747,assembly=hg19>\n\
##contig=<ID=chr11,length=135006516,assembly=hg19>\n\
##contig=<ID=chr12,length=133851895,assembly=hg19>\n\
##contig=<ID=chr13,length=115169878,assembly=hg19>\n\
##contig=<ID=chr14,length=107349540,assembly=hg19>\n\
##contig=<ID=chr15,length=102531392,assembly=hg19>\n\
##contig=<ID=chr16,length=90354753,assembly=hg19>\n\
##contig=<ID=chr17,length=81195210,assembly=hg19>\n\
##contig=<ID=chr18,length=78077248,assembly=hg19>\n\
##contig=<ID=chr19,length=59128983,assembly=hg19>\n\
##contig=<ID=chr20,length=63025520,assembly=hg19>\n\
##contig=<ID=chr21,length=48129895,assembly=hg19>\n\
##contig=<ID=chr22,length=51304566,assembly=hg19>\n\
##contig=<ID=chrX,length=155270560,assembly=hg19>\n\
##contig=<ID=chrY,length=59373566,assembly=hg19>\n'


# sort all sampleIDs, just in case
all_names.sort()

# open conversion chart generated by hgvsC_to_vcf-like.py
# has hgvsC and resulting vcf-like output (chr, pos, ref, alt)

# rolling through each line in vcf now to reformat and print to file
file = open('gemini.casey.vcf','w')
for line in vcf:
	# reprint header
	if line[0:2] == '##':
		file.write(line)
	# first add new header lines, then print amended column names (FORMAT + sample names)
	elif line[0:2] == '#C':
		file.write(new_header)
		line = line.split()
		line.append("FORMAT")
		line = line + all_names
		file.write('\t'.join(line))
		file.write('\n')
	# now out of header. Now hacking vcf with sample genotypes and fake AD:DP:PL (see below)
	else:
		line = line.split()
		# hacking conversion from ensembl chr to ucsc chr
		line[0] = "chr" + line[0]
		# Faking AD:DP:PL (see below)
		line.append("GT:AD:DP:PL")
		# pull genotypes at the position by running HGVS through dict
		genotypes = hgvs_dict[line[2]]
		for sample in all_names:
			if sample not in genotypes:
				# fake AD:DP:PL numbers to appease Gemini
				line.append("./.:666:666:666")
			else:
				# casey uses some version of "heterozygote" to mark hets for an allele
				if 'het' in hgvs_dict[line[2]][sample].lower():
					line.append("0/1:666,666:666:666")
				# similar, homozygous for hom for alt allele
				elif 'hom' in hgvs_dict[line[2]][sample].lower():
					line.append("1/1:666:666:666")
				# possible something weird might slip in. Kill script and alert user.
				else:
					print('Error, what is: ', hgvs_dict[line[2]][sample].lower())
					print('\n\nKILL KILL KILL\n\n')
					break
		file.write('\t'.join(line))
		file.write('\n')

file.close()			

##########################################################################################
# gemini db creation
##########################################################################################
# attaching current date and time to db for rudimentary version control
db_name = "casey.gemini." + datetime.now().strftime('%Y-%m-%d__%H:%M:%S') + ".db"
gemini_query = "gemini load --cores 4 -t VEP -v gemini.casey.vcf -p casey.ped " + db_name

os.system(gemini_query)
		




