from os.path import join


configfile: '/home/mcgaugheyd/git/casey_to_gemini/config.yaml'


def find_excel_sheet_match(wildcards):
	import glob
	sample = str(wildcards).split('_')[0]
	excel_file = glob.glob('excel_report/' + sample + '*')[0]
	return(excel_file)

(ALL_SAMPLES, END) = glob_wildcards(join('vcf/', '{sample}.vc{end}'))
(NEW_SAMPLES, ) = glob_wildcards(join('vcf/', '{sample}.vcf'))
(EXISTING_SAMPLES, ) = glob_wildcards(join('vcf/', '{sample}.vcf.gz'))

localrules: move_existing_vcf_to_temp, MVL_excel_sheet_query


rule all:
	input:
		config['gemini_db_name'],
		expand('MGOG_reports/{sample}.report.html', sample=NEW_SAMPLES)


rule vt_bgzip_and_tabix_MVL_vcf:
	input:
		'vcf/{sample}.vcf'
	output:
		vcf = temp('temp/{sample}.vcf.gz'),
		index = temp('temp/{sample}.vcf.gz.tbi')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		# hack for now, replace M with MT
		cat {input} \
			| sed 's/^M/MT/g' - \
			| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
			| grep -v "##sgmutationstatistics" \
			| ~/git/vt/./vt decompose -s - \
			| ~/git/vt/./vt normalize -n -r {config[ref_genome]} - \
			| bgzip -c > {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""

rule merge_all_vcfs:
	input:
		vcf_new = expand('temp/{sample}.vcf.gz', sample=NEW_SAMPLES),
		vcf_existing = expand('vcf/{sample}.vcf.gz', sample=EXISTING_SAMPLES),
		index = expand('temp/{sample}.vcf.gz.tbi', sample=NEW_SAMPLES),
		index_existing = expand('vcf/{sample}.vcf.gz.tbi', sample=EXISTING_SAMPLES)
	output:
		temp('temp/ALL_SAMPLES.vcf')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		bcftools merge {input.vcf_new} {input.vcf_existing} > {output}
		"""

rule sort_bgzip_and_tabix_big_vcf:
	input:
		'temp/ALL_SAMPLES.vcf'
	output:
		vcf = temp('temp/ALL_SAMPLES.SORTED.vcf.gz'),
		index = temp('temp/ALL_SAMPLES.SORTED.vcf.gz.tbi')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		mkdir -p /scratch/$SLURM_JOB_ID/
		grep '^#' {input} > /scratch/$SLURM_JOB_ID/TEMPTEMP
		grep -v '^#' {input} | sort -k1,1 -k2,2n >> /scratch/$SLURM_JOB_ID/TEMPTEMP 
		bgzip /scratch/$SLURM_JOB_ID/TEMPTEMP -c > {output.vcf}
		tabix -f -p vcf {output.vcf}
		"""

rule make_fake_ped:
	input:
		vcf = 'temp/ALL_SAMPLES.SORTED.vcf.gz',
		index = 'temp/ALL_SAMPLES.SORTED.vcf.gz.tbi'
	output:
		temp('temp/FAKE.ped')
	shell:
		"""
		module load {config[samtools_version]}
		samples=`bcftools query -l {input.vcf}`
		for i in $samples; do 
			echo $i $i 0 0 0 2; 
		done > {output}
		"""

rule vt_left_align_master_vcf:
	input:
		vcf = 'temp/ALL_SAMPLES.SORTED.vcf.gz',
		index = 'temp/ALL_SAMPLES.SORTED.vcf.gz.tbi'
	output:
		temp('temp/ALL_SAMPLES.SORTED.VT.vcf.gz')
	shell:
		"""
		cat {input.vcf} \
			| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
			| ~/git/vt/./vt decompose -s - \
			| ~/git/vt/./vt normalize -n -r {config[ref_genome]} - \
			> {output}
		"""

rule VEP_annotate:
	input:
		'temp/ALL_SAMPLES.SORTED.VT.vcf.gz'
	output:
		vcf = temp('temp/ALL_SAMPLES.SORTED.VT.VEP.vcf.gz'),
		index = temp('temp/ALL_SAMPLES.SORTED.VT.VEP.vcf.gz.tbi')
	threads: 16
	shell:
		"""
		module load {config[VEP_version]}
		vep -i {input} --offline \
			--cache --dir_cache $VEPCACHEDIR \
			--fasta $VEPCACHEDIR/GRCh37.fa --species human --assembly GRCh37  \
			--format vcf \
			--output_file {output.vcf} \
			--plugin Grantham \
			--plugin MaxEntScan,/data/OGVFB/resources/MaxEntScan \
			--plugin CADD,/fdb/CADD/1.3/prescored/whole_genome_SNVs.tsv.gz,/fdb/CADD/1.3/prescored/InDels.tsv.gz \
			--canonical \
			--ccds \
			--total_length \
			--hgvs \
			--sift b \
			--polyphen b \
			--symbol \
			--numbers \
			--biotype \
			--total_length \
			--pubmed \
			--domains \
			--gene_phenotype \
			--pick \
			--pick_order canonical, tsl, biotype, ccds, length \
			--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,DOMAINS,CLIN_SIG,Grantham,MaxEntScan,HGVSc,HGVSp,PUBMED,Phenotypes,CADD_RAW,CADD_PHRED \
			--vcf --compress_output bgzip --force_overwrite --fork {threads}
		# tabix
		tabix -f -p vcf {output.vcf}
		"""

rule vcfanno_annotate:
	input:
		vcf = 'temp/ALL_SAMPLES.SORTED.VT.VEP.vcf.gz',
		index = 'temp/ALL_SAMPLES.SORTED.VT.VEP.vcf.gz'
	output:
		vcf = temp('temp/ALL_SAMPLES.SORTED.VT.VEP.VCFANNO.vcf.gz'),
		index = temp('temp/ALL_SAMPLES.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi')
	threads: 16
	shell:
		"""
		module load {config[vcfanno_version]}
		vcfanno -p {threads} -lua {config[vcfanno_lua]} {config[vcfanno_conf]} {input.vcf} | bgzip > {output.vcf} 
		tabix -f -p vcf {output.vcf}
		"""

rule make_gemini_db:
	input:
		vcf = 'temp/ALL_SAMPLES.SORTED.VT.VEP.VCFANNO.vcf.gz',
		index = 'temp/ALL_SAMPLES.SORTED.VT.VEP.VCFANNO.vcf.gz.tbi',
		ped = 'temp/FAKE.ped'
	output:
		config['gemini_db_name']
	shell:
		"""
		module load {config[vcf2db_version]}
		vcf2db.py {input.vcf} {input.ped} {output}
		"""

rule query_gemini:
	input:
		db = config['gemini_db_name'],
		input_vcf = 'vcf/{sample}.vcf'
	output:
		temp('temp/{sample}/annotated_variants')
	shell:
		"""
		export REF_CACHE=/scratch/$SLURM_JOB_ID/
		module load {config[samtools_version]}
		sample=`bcftools query -l {input.input_vcf} | sed 's/-/_/g'`
		module load {config[gemini_version]}
		gemini query --show-samples --header -q "select chrom, start, end, ref, alt, vcf_id, hgvsc, hgvsp, type, gene, impact, impact_severity, rs_ids, pfam_domain, clinvar_diseases, clinvar_ID, clinvar_pathogenic, clinvar_sig, hgmd_overlap, num_het, num_hom_alt, max_aaf_all, af_exac_afr, af_exac_all, af_exac_amr, af_exac_eas, af_exac_nfe, af_exac_oth, af_exac_sas, an_exac_all, exac_num_het, exac_num_hom_alt, mis_z, pli, gerp_elements, polyphen_pred, polyphen_score, sift_pred, sift_score, cadd_phred, gt_depths.$sample, gts.$sample, gt_quals.$sample FROM variants" --gt-filter "gt_types.$sample == HET or gt_types.$sample == HOM_ALT" {input.db} > {output}  
		"""

rule MVL_excel_sheet_query:
	input:
		find_excel_sheet_match
	output:
		temp('temp/{sample}/mvl_sheet_info')
	shell:
		"""
		/usr/local/Anaconda/envs/py3.4/bin/python3 /home/mcgaugheyd/git/casey_to_gemini/src/casey_xlsx_to_hgvs.py {input} > {output}
		"""

rule write_report:
	input:
		annotated_variants = 'temp/{sample}/annotated_variants',
		mvl_excel_info = 'temp/{sample}/mvl_sheet_info'
	output:
		'MGOG_reports/{sample}.report.html'
	shell:
		"""
		module load {config[R_version]}
		module load {config[pandoc_version]}
		cp /home/mcgaugheyd/git/casey_to_gemini/src/write_report.Rmd temp/{wildcards.sample}/
		cp /home/mcgaugheyd/git/casey_to_gemini/src/write_report.R temp/{wildcards.sample}/
		cp /home/mcgaugheyd/git/casey_to_gemini/src/datatable_fix.css temp/{wildcards.sample}/
		Rscript temp/{wildcards.sample}/write_report.R \
			temp/{wildcards.sample}/write_report.Rmd \
			annotated_variants \
			mvl_sheet_info \
			{wildcards.sample} \
			report.html
		mv temp/{wildcards.sample}/report.html {output}
		rm temp/{wildcards.sample}/write_report.Rmd
		rm temp/{wildcards.sample}/write_report.R
		rm temp/{wildcards.sample}/datatable_fix.css
		"""
