#!/bin/bash

# show phenotype info available
gemini query --header -q "SELECT * FROM samples" casey.gemini.2016-02-18__16\:16\:21.db


# shows all info you can search on
gemini db_info casey.gemini.2016-02-18__16\:16\:21.db | less -S


#Adds up number of variants (hom_ref, het, hom_alt, across each person)
gemini stats --summarize "select * from variants" casey.gemini.2016-02-18__16\:16\:21.db | csvlook -t


#similar, but summaries only predicted clinvar pathogenic variants
gemini stats --summarize "select * from variants WHERE clinvar_sig LIKE '%path%'" casey.gemini.2016-02-18__16\:16\:21.db | csvlook -t


#or predicted HIGH biological impact
gemini stats --summarize "select * from variants WHERE impact_severity='HIGH'" casey.gemini.2016-02-18__16\:16\:21.db | csvlook -t


#patient B_B comes up with a few, so let's see what variants are predicted pathogenic for B_B
pat_name='patient_name'

gemini query --header \
-q "SELECT gene, chrom, start, ref, alt, name, Disease, (gts).(name==pat_name),vep_hgvsc, vep_hgvsp FROM variants,samples \
WHERE name=pat_name AND impact_severity='HIGH'" \
casey.gemini.2016-02-18__16\:16\:21.db| csvlook -t | less -S


# Summarize people with a particular variant and stratify by a phenotype (Disease)
gemini query --show-samples --carrier-summary-by-phenotype Disease --header -q "select gene, chrom, start, ref, alt, vep_hgvsc, vep_hgvsp, (gts).(*) from variants" casey.gemini.2016-02-18__16\:16\:21.db | csvsort -t -c 1 | csvlook | less -S


# show variants with achromotopsia subjects with non-ref alleles ("!=HOM_REF", != means "not equals")
# can also change to just homozygous alt ("==HOM_ALT")
# or just het ("==HET")
# only showing genotypes for people with Achromotopsia ("(gts).(Disease=='Achromotopsia')")
gemini query --header \
-q "SELECT gene, chrom, start, ref, alt, name, Disease, (gts).(Disease=='Achromotopsia'),vep_hgvsc, vep_hgvsp \
FROM variants,samples WHERE Disease='Achromotopsia'" \
--gt-filter "(gt_types).(Disease=='Achromotopsia').(!=HOM_REF).(any)" \
casey.gemini.2016-02-18__16\:16\:21.db | csvlook -t | less -S


# same, but where 2 or more patients have the same genotype ("count>=2")
gemini query --header \
-q "SELECT gene, chrom, start, ref, alt, name, Disease, (gts).(Disease=='Achromotopsia'),vep_hgvsc, vep_hgvsp \
FROM variants,samples WHERE Disease='Achromotopsia'" --gt-filter "(gt_types).(Disease=='Achromotopsia').(!=HOM_REF).(count>=2)" \
casey.gemini.2016-02-18__16\:16\:21.db | csvlook -t | less -S
