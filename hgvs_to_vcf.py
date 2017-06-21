#!/usr/local/bin/python2

import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import sys
import subprocess

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
		# fix output for deletion, since you need the base 5' to the HGVS given
		# vcf for deletion takes the format chr, pos (1 before the deleletion), 
			# ref (including the 1 base pair before deletion), alt (just the 1 base pair before deletion)
		if alt is None and 'del' in hgvs_c:
			if reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos+2).upper() == ref:
				pos = str(int(pos)-1)
				ref = reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos).upper() + ref
				alt = reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos).upper() 
			else:
				alt = '.'
		output = [str(chr), str(pos), '.', str(ref), str(alt)]
		return output
	except Exception as e:
		output = ["ERROR: ", str(e)]
		return(output)

def reference(fasta, chrom, start, stop):
	""" Take path to fasta file, then use samtools faidx to gragb reference base in region """
	region = chrom + ':' + str(start) + '-' + str(stop)
	# fasta = '/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta'
	call = 'samtools faidx ' + fasta + ' ' + region
	output = subprocess.check_output(call, shell=True)
	return(output.split('\n')[1])

if __name__ == '__main__':
	vcf_like = converter(sys.argv[1])
	print('\t'.join(vcf_like))
