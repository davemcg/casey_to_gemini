#!/usr/local/bin/python3

import pandas
import math
import os
import subprocess
import re
from datetime import datetime
from collections import defaultdict
import sys

##########################################################################################
#Reads in xlsx file, extract hgvs, and writes ped file
##########################################################################################

# Reads all files in directory
all_files = os.listdir()
# only input patient xlsx files and builds a ped file (doesn't do gender)
xlsx = []
ped = []
ped.append('#Family SampleID PaternalID MaternalID Gender Phenotype ClinicalID')

#build list of all samples for use through script
all_names = []
for i in all_files:
	# ID file names that start with a digit and end in xlsx
	if re.search('^\d+.*xlsx$',str(i)):
		xlsx.append(i)
		i = i.replace('_',' ')
		clinicalID = i.split()[0]
		# extract first and last name from file name
		fname = i.split()[1]		
		lname = i.split()[2]
		name = fname + '_' + lname
		all_names.append(name)
		# and pull disease from file name
		# disease = i.split()[3] #NOPE TOO SKETCHY. FILENAMING INCONSISTENT.
		# clean up a formatting thing where some diseases are appended with "_"
		#if '_' in disease:
		#	disease = disease.split('_')[0]
		ped_line = name + ' ' + name + ' 0 ' + '0 ' + '0 ' + '1 ' + clinicalID
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

# loops through xlsx list, opens xlsx file to extract hgvs, build three dictionaries: 
# 1: position_zygosity_dict, which has position (hgvsCoding) as key and zygosity as value
# 2: position_panel_dict, position (hgvsCoding) as key and panel used as value
# 3: sample_panel_dict, sample (subject/person) as key and panel used as value
panels = []
position_zygosity_dict = defaultdict(dict)
position_panel_dict = defaultdict(dict)
sample_panel_dict = defaultdict(dict)
for i in xlsx:
#	print('Rolling through ' + i + ' now')
	## check to find how many rows to skip and grab panel info 
	check = pandas.read_excel(i,sheetname=0)
	# pull first column and convert to list, then check for row with #ID (header row)
	first_column = list(check[check.columns[0]])
	index = first_column.index('#ID') + 1
	# the index happens to be position of the Sample and panel info cell
	panel_index = [i for i,val in enumerate(first_column) if str(val).startswith("##Variant")][0]
	panel = first_column[panel_index]
	panel = re.search("(?<=\d\d\d[-|\s+]).*\'", panel).group(0)
	panel = panel.replace('\'','')
	panel = panel.replace('-','_')
	panel = panel.replace(' ','_')
	print('Going through ' + i + ' | ' + panel)
	# read first xlsx sheet into panda data structure 
	data = pandas.read_excel(i, sheetname=0, skiprows=index)
	#just keep select info, depending on what kind of sheet I get
	if 'Zygosity' not in data:
		print('Zygosity Information missing from ' + i)
		# move xlsx to problem_xlsx_files
		if not os.path.isdir("problem_xlsx_files"):
			mkdir_call = "mkdir problem_xlsx_files"
			subprocess.check_call(mkdir_call,shell=True)
		mv_call = 'mv ' + '\"' + i + '\"' + ' problem_xlsx_files/'
		subprocess.check_call(mv_call,shell=True)
		continue
	if 'Transcript' in data: # early sheets have transcript and hgvs in separate columns
		data = data[['Transcript','HGVSCoding','Zygosity','TimesObservedPerPanel']]
		# merge back together
		data['HGVSCoding'] = str(data['Transcript'].map(str)) + ':' + str(data['HGVSCoding'])
	else:
		data = data[['HGVSCoding','Zygosity','TimesObservedPerPanelGroup','SamplesPerPanelGroup']]
	# drop blanks (in any column)
	data = data.dropna()
	# drop anything that doesn't have a HGVS in it
	cond = ~data['HGVSCoding'].str.contains('c.')
	data = data.drop(data[cond].index.values)
	# build dictionaries
	fname = i.split()[1]		
	lname = i.split()[2]
	name = fname + "_" + lname
	for index,row in data.iterrows():
		key = row['HGVSCoding']
		if 'M_' not in key[0:3]:
			continue
		if ':' not in key:
			print(key, "something wrong???")
			sys.exit(0)
		key1, key2 = key.split(':')
		key1 = key1.split('.')[0]
		# warn about the occasional missing (MT?)
		if key1 == ' ' or key1 == '':
			print(i, key1, key2, "Missing gene!!!!!!\nNot being entered into db!!!!!\n\n")
			continue
		if key1 == 'NM_000114':
			# NM_000114 is gone now. Replacing with below
			key1 = 'NM_207034'
		if key1 == 'XM_005273027':
			key1 = 'NM_001301365'
		key = key1 + ':' + key2
		panels.append(panel)
		position_zygosity_dict[key][name] = row['Zygosity'] + ";;" + str(row['TimesObservedPerPanelGroup']) + ";;" + str(row['SamplesPerPanelGroup'])
		if key not in position_panel_dict:	
			position_panel_dict[key] = panel
		if key in position_panel_dict:
			value = list(position_panel_dict[key])
			value.append(panel)
		if name not in sample_panel_dict:
			sample_panel_dict[name] = panel
		if name in sample_panel_dict:
			value = list(sample_panel_dict[name])
			value.append(panel)
			value = set(value)


print("\n\nHere are the panels used across the files: ", set(position_panel_dict.values()), "\n\n")
# print hgvs to file for counsyl hvgs (hgvsC_to_vcf-like.py)
file=open("hgvs_input.txt",'w')
for k,v in position_zygosity_dict.items():
	file.write(k)
	file.write('\n')

file.close()

##########################################################################################
# Convert to VEP friendly format and re-do position_zygosity_dict with chr_pos_ref_alt as key
##########################################################################################
os.system('~/git/casey_to_gemini/./hgvsC_to_vcf-like.py hgvs_input.txt')
converter = open('converter.txt')
for line in converter:
	line = line.split()
	# update position_zygosity_dict
	values = position_zygosity_dict[line[0]]
	new_key = 'chr' + line[1]
	position_zygosity_dict[new_key] = values
	# update position_panel_dict
	panel_values = position_panel_dict[line[0]]
	position_panel_dict[new_key] = panel_values

# Change 

##########################################################################################
# Run VEP. 
##########################################################################################

# will output vcf from the hgvs
vep_query = 'perl /Applications/variant_effect_predictor/variant_effect_predictor.pl ' + \
		'-i vep_hgvs_input.txt ' + \
		'--vcf -o ' + 'vep_output.vcf ' + \
		'--assembly GRCh37 ' + \
		'--offline ' + \
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

print("Running VEP. Expect to wait ~10 minutes")
subprocess.check_call(vep_query,shell=True)
print("VEP query done!")

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
##INFO=<ID=caseyTOPPG,Number=R,Type=Integer,Description="Casey TimesObservedPerPanelGroup">\n\
##INFO=<ID=caseySPPG,Number=R,Type=Integer,Description="Casey SamplesPerPanelGroup">\n\
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
	#
	# now out of header. Now hacking vcf with sample genotypes and fake AD:DP:PL (see below)
	#
	else:
		line = line.split()
		# hacking conversion from ensembl chr to ucsc chr
		line[0] = "chr" + line[0]
		# Faking AD:DP:PL (see below)
		line.append("GT:AD:DP:PL")
		# build vcf_key from vcf info
		vcf_key = line[0] + '_' + line[1] + '_' + line[3] + '_' + line[4]

		# pull genotypes at the position by running HGVS through dict
		sample_genotypes = position_zygosity_dict[vcf_key]
		for sample in all_names:
			if sample not in sample_genotypes:
				# begin the craziness. Need to call a db to check:
				# 1. what panel is used on this sample/subject?
				# 2. is this position covered on the panel?
				# if 2 is true, then mark position as 0/0 (hom_ref)
				# otherwise mark at unknown (./.)
				if sample_panel_dict[sample] == position_panel_dict[vcf_key]:
					line.append("0/0:666:666:666")
				else:
					# fake AD:DP:PL numbers to appease Gemini
					line.append("./.:666:666:666")
			else:
				line[7] = "caseyTOPPG=" + position_zygosity_dict[vcf_key][sample].split(';;')[1] + '|' + "caseySPPG=" + position_zygosity_dict[vcf_key][sample].split(';;')[2]
				# casey uses some version of "heterozygote" to mark hets for an allele
				if 'het' in position_zygosity_dict[vcf_key][sample].split(';;')[0].lower():
					line.append("0/1:666:666:666")
				# similar, homozygous for hom for alt allele
				elif 'hom' in position_zygosity_dict[vcf_key][sample].split(';;')[0].lower():
					line.append("1/1:666:666:666")
				# possible something weird might slip in. Alert user.
				else:
					line.append('./.:000:000:000')
					print('Error, what is: ', sample, position_zygosity_dict[vcf_key][sample].split(';;')[0].lower())
					print('Genotype unknown: ./.')
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
		




