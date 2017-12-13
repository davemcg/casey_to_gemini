#!/usr/local/bin/python2

import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import sys
import subprocess
import argparse

# connect to mapping database
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
vm37 = hgvs.assemblymapper.AssemblyMapper(
    hdp, assembly_name='GRCh37', alt_aln_method='splign')
vr = hgvs.validator.Validator(hdp=hdp)

# unit test inputs
hgvs_c = 'NM_018474.4:c.226C>T'


def converter(hgvs_c):
	""" Take HGVS (coding) input and return chr, GRCh37 position, reference allele, alternative allele """
	
	# parse into hgvs structure
	try:
		var_c = hp.parse_hgvs_variant(hgvs_c)
		vr.validate(var_c)
	except Exception as e:
		output = ["ERROR: ", str(e), hgvs_c]
		return(output)

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
			if reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos+len(ref)-1).upper() == ref:
				pos = str(int(pos)-1)
				ref = reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos).upper() + ref
				alt = reference('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta', 'chr' + chr, pos, pos).upper() 
			else:
				raise ValueError('Failed to convert deletion {hgvs}'.format(hgvs=repr(hgvs_c)))
		output = [str(chr), str(pos), str(hgvs_c), str(ref), str(alt)]
		return output
	except Exception as e:
		output = ["ERROR: ", str(e), hgvs_c]
		return(output)

def reference(fasta, chrom, start, stop):
	""" Take path to fasta file, then use samtools faidx to grab reference base in region """
	region = chrom + ':' + str(start) + '-' + str(stop)
	# fasta = '/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta'
	call = 'samtools faidx ' + fasta + ' ' + region
	output = subprocess.check_output(call, shell=True)
	return(output.split('\n')[1])

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description= 'Takes in c. HGVS and converts to a vcf like format')
	parser.add_argument('--single','-s', help= 'Single HGVS input')
	parser.add_argument('--file', '-f', help= 'New line separated file of HGVS')
	parser.add_argument('--comma', '-c', help= 'Comma separated c. HGVS')
	args = parser.parse_args()

	if args.single:
		vcf_like = converter(args.single)
		print('\t'.join(vcf_like))
	if args.file:
		hgvs_info = open(args.file,'r')
		for line in hgvs_info:
			vcf_like = converter(line[:-1])
			print('\t'.join(vcf_like))
	if args.comma:
		hgvs_info = args.comma.split(',')
		for line in hgvs_info:
			vcf_like = converter(line)
			print('\t'.join(vcf_like))
