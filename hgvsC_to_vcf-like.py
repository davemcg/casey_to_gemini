#!/usr/local/bin/python2

# counsyl hgvs

import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
import sys
input = sys.argv[1]

# Read genome sequence using pygr.
from pygr.seqdb import SequenceFileDB
genome = SequenceFileDB('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta')
# Read RefSeq transcripts into a python dict.
with open('/Users/mcgaugheyd/GenomicData/hg19.refGene') as infile:
	transcripts = hgvs_utils.read_transcripts(infile)
# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
	return transcripts.get(name)


c = open('converter.txt','w')
f = open('vep_hgvs_input.txt','w')
for line in open(input):
	chrom, offset, ref, alt = hgvs.parse_hgvs_name(line, genome, get_transcript=get_transcript)
	chrom = chrom[3:]
	out = str(chrom) + ' ' + str(offset) + ' . ' + ref + ' ' + alt
	f.write(out)
	out2 = line[:-1] + '\t' + str(chrom) + '_' + str(offset) + '_' + ref + '_' + alt
	c.write(out2)
f.close()
c.close()
