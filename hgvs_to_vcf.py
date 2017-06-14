#!/usr/local/bin/python2

import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import sys

# connect to mapping database
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
vm37 = hgvs.assemblymapper.AssemblyMapper(
    hdp, assembly_name='GRCh37', alt_aln_method='splign')

# unit test inputs
hgvs_c = 'NM_018474.4:c.226C>T'


def converter(hgvs_c):
	""" Take HGVS (coding) input and return chr, GRCh37 position, reference allele, alternative allele """
	
	# parse into hgvs structure
	var_c = hp.parse_hgvs_variant(hgvs_c)
	# translate into genomic coordinates
	try:
		var_g = vm37.c_to_g(var_c)
		# chr
		chr = str(int(var_g.ac.split('.')[0].split('_')[1]))
		if chr == '23':
			chr = 'X'
		if chr == '24':
			chr = 'Y'
		# position
		pos = var_g.posedit.pos.start.base
		# ref
		ref = var_g.posedit.edit.ref
		# alt
		alt = var_g.posedit.edit.alt

		output = [str(chr), str(pos), '.', str(ref), str(alt)]
		return output
	except Exception as e:
		output = ["ERROR: ", str(e)]
		return(output)

if __name__ == '__main__':
	vcf_like = converter(sys.argv[1])
	print('\t'.join(vcf_like))
