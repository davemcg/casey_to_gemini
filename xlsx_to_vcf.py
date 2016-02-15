import pandas
import math
import sys

"""
Read in xlsx file, extract hgvs, and ped info
"""

x = pandas.read_excel("/Users/davidmcgaughey/Desktop/casey/15-3081 Victoria Sang Stargardt v3.2_Report-sign out.xlsx", sheetname=0,skiprows=4)

#list(x.columns.values)
#['#ID', 'Gene', 'Chromosome', 'ChromosomePosition', 'ExonNumber', 'HGVSGenomic', 'HGVSCoding', 'HGVSProtein', 'Zygosity', 'Coverage', 'Rs', 'CAF[RA]', 'ExAC_AF', 'Pathogenicity', 'Polyphen2_HDIV_pred', 'Polyphen2_HDIV_score', 'VariantComment', 'TimesObservedPerPanel', 'SamplesPerPanel', 'TimesObservedPerPanelGroup', 'SamplesPerPanelGroup', 'Panel', 'Unnamed: 22', 'Unnamed: 23', 'Unnamed: 24', 'Unnamed: 25']

#x['HGVSGenomic']

#y = x[['Chromosome','HGVSGenomic']]

x['HGVS'] = x["Chromosome"].map(str) + ":" + x["HGVSGenomic"]

out = x['HGVS'].tolist()

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
		




