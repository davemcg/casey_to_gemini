import pandas

x = pandas.read_excel("/Users/davidmcgaughey/Desktop/casey/15-3081 Victoria Sang Stargardt v3.2_Report-sign out.xlsx", sheetname=0,skiprows=4)

list(x.columns.values)
['#ID', 'Gene', 'Chromosome', 'ChromosomePosition', 'ExonNumber', 'HGVSGenomic', 'HGVSCoding', 'HGVSProtein', 'Zygosity', 'Coverage', 'Rs', 'CAF[RA]', 'ExAC_AF', 'Pathogenicity', 'Polyphen2_HDIV_pred', 'Polyphen2_HDIV_score', 'VariantComment', 'TimesObservedPerPanel', 'SamplesPerPanel', 'TimesObservedPerPanelGroup', 'SamplesPerPanelGroup', 'Panel', 'Unnamed: 22', 'Unnamed: 23', 'Unnamed: 24', 'Unnamed: 25']

x['HGVSGenomic']

y = x[['Chromosome','HGVSGenomic']]

y['HGVS'] = y["Chromosome"].map(str) + ":" + y["HGVSGenomic"]