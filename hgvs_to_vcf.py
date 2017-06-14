#!/usr/bin/env python2

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
