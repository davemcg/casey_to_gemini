"""
Invitae HGVS
hgvs.readthedocs.org

pros:
	- in active development
	- error/sanity checks
cons:
	- can't do VCF style yet
	- relies on external dataset (UTA)
		- but perhaps could be installed locally?
"""

import hgvs.parser
import hgvs.variantmapper
import hgvs.dataproviders.uta

# our variant
hgvs_c = 'NM_000350.2:c.4774-17_4774-16delGT'

# connect to mapping database
hdp = hgvs.dataproviders.uta.connect()
evm = hgvs.variantmapper.EasyVariantMapper(hdp,primary_assembly='GRCh37', alt_aln_method='splign')

# parse into hgvs structure
hp = hgvs.parser.Parser()
var_c = hp.parse_hgvs_variant(hgvs_c)
# translate into genomic coordinates
var_g = evm.c_to_g(var_c)


# start
var_g.posedit.pos.start.base
# ref
var_g.posedit.edit.ref
# alt
var_g.posedit.edit.alt

##########################################################################################
"""
Counsynl hgvs
www.github.com/counsyl/hgvs

pros: 
	- works now
	- uses local files
negs:
	- not in active development?
	- refgene file doesn't do transcripts (i.e. NM_000352.3 is just called NM_000352. A potential issue?)
	- doesn't error check (won't return error if ref not matched)
	
"""

import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils
from pygr.seqdb import SequenceFileDB

# Read genome sequence using pygr.
genome = SequenceFileDB('/Users/mcgaugheyd/GenomicData/ucsc.hg19.fasta')

# Read RefSeq transcripts into a python dict.
with open('/Users/mcgaugheyd/GenomicData/hg19.refGene') as infile:
    transcripts = hgvs_utils.read_transcripts(infile)

# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
	return transcripts.get(name)

# Parse the HGVS name into genomic coordinates and alleles.
chrom, offset, ref, alt = hgvs.parse_hgvs_name('NM_000352:c.215A>G', genome, get_transcript=get_transcript)
# Returns variant in VCF style: ('chr11', 17496508, 'T', 'C')
# Notice that since the transcript is on the negative strand, the alleles
# are reverse complemented during conversion.

# Format an HGVS name.
chrom, offset, ref, alt = ('chr11', 17496508, 'T', 'C')
transcript = get_transcript('NM_000352.3')
hgvs_name = hgvs.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript)
# Returns 'NM_000352.3(ABCC8):c.215A>G'