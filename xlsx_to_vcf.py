import pandas
import math
import os

"""
Read in xlsx file, extract hgvs, and build ped
"""

# Will design right now to read all xlsx files in directory
all_files = os.listdir()
# only read in patient xlsx files and builds a ped file (doesn't do gender)
xlsx = []
ped = []
ped.append('#Family SampleID PaternalID MaternalID Gender Phenotype Disease ClinicalID')
for i in files:
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
	data = pandas.read_excel(i, sheetname=0, skiprows=4)
	data['HGVS'] = data["Chromosome"].map(str) + ":" + data["HGVSGenomic"]
	out = data['HGVS'].tolist()
	# filter out non hgvs rows
	out2 = [] 
	for object in out:
		if re.search('^\d+', str(object)):
			out2.append(object)
		
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
		




