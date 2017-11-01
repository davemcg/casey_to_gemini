#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from openpyxl import load_workbook
import openpyxl
import argparse
import subprocess
import datetime
import time

parser = argparse.ArgumentParser(description= 'Run on eyeMac (local computer). Takes in xlsx files from John Chiang / MVLGenomics Gene Panels and does the following:\n 1. Identifies hidden rows and skips these \n 2. Extracts hgvs \n\n Usage: casey_xlsx_to_hgvs.py YOURFILE.xlsx > output.tsv')
parser.add_argument('xlsx_file', help = 'John Chiang / OSU / MVL / Casey patient report xlsx file')

args = parser.parse_args()
print('Begin processing of ' + args.xlsx_file.split(' ')[0] + ' at ' + str(datetime.datetime.now()))
# load workbook and read sheet names with openpyxl
wb = load_workbook(args.xlsx_file)
ws = wb[wb.get_sheet_names()[0]]
# identify header row, HGVS coding column, and zygosity column
for row in ws.iter_rows(min_row=2, max_row=10):
    for cell in row:
        if cell.value == 'HGVSCoding':
            header_row = cell.row
            HGVS_column = cell.column
        if cell.value == 'Zygosity':
            Zygosity_column = cell.column
        if cell.value == 'Panel':
            Panel_column = cell.column
        if cell.value == 'Transcript':
            transcript_column = cell.column
# identify last row (John sometimes puts crap at the bottom)
for row in range(header_row+1, ws.max_row):
    if ws[HGVS_column + str(row)].value is None:
        last_row = row
        break

# get panel type
try:
    for row in range(header_row+1, ws.max_row):
        if ws[Panel_column + str(row)].value is not None:
            panel_type = ws[Panel_column + str(row)].value.replace(' ','_').replace('(','_').replace(')','_')
except:
    panel_type = '_'.join(args.xlsx_file.split(' ')[3:])
    panel_type = panel_type.replace('.xlsx','').replace('(','_').replace(')','_')

# fill in value for transcript_column if missing (only OLD 2015 panels have this field)
try:
    transcript_column
except:
    transcript_column = 'D'

# get HGVS and zygosity
hgvs_zygosity = {}
hgvs_status = {}
hgvs_vars = []
for row in range(header_row+1, last_row):
    cell_name_hgvs = "{}{}".format(HGVS_column, row)
    cell_name_zyg = "{}{}".format(Zygosity_column, row)
    cell_name_tx = "{}{}".format(transcript_column, row)
    if ws.row_dimensions[row].hidden is False: # not hidden
        if ws[cell_name_hgvs].value is not None and ws[cell_name_tx].value is not None and (ws[cell_name_hgvs].value[0:2] == 'NM' or ws[cell_name_tx].value[0:2] == 'NM'):
            if 'NM' not in ws[cell_name_hgvs].value:
                 hgvs_zygosity[ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value] = ws[cell_name_zyg].value
                 hgvs_status[ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value] = 'Primary'
                 hgvs_vars.append(ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value)
            else:
                hgvs_zygosity[ws[cell_name_hgvs].value] = ws[cell_name_zyg].value
                hgvs_status[ws[cell_name_hgvs].value] = 'Primary'
                hgvs_vars.append(ws[cell_name_hgvs].value)
    else: # is hidden
        if ws[cell_name_hgvs].value is not None and ws[cell_name_tx].value is not None and (ws[cell_name_hgvs].value[0:2] == 'NM' or ws[cell_name_tx].value[0:2] == 'NM'):
            if 'NM' not in ws[cell_name_hgvs].value:
                hgvs_zygosity[ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value] = ws[cell_name_zyg].value
                hgvs_status[ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value] = 'Secondary'
                hgvs_vars.append(ws[cell_name_tx].value + ':' + ws[cell_name_hgvs].value)
            else:
                hgvs_zygosity[ws[cell_name_hgvs].value] = ws[cell_name_zyg].value
                hgvs_status[ws[cell_name_hgvs].value] = 'Secondary'
                hgvs_vars.append(ws[cell_name_hgvs].value)
# first attempt local conversion
print('Local Conversion Begun for ' + args.xlsx_file.split(' ')[0])
local_conversion_results = subprocess.check_output(['/Users/mcgaugheyd/git/casey_to_gemini/src/hgvs_to_vcf.py','--comma', ','.join(hgvs_vars)]).decode('utf-8')
print('VEP Conversion Begun for ' + args.xlsx_file.split(' ')[0])
# then take the failures and run against VEP
hgvs_file_name =  args.xlsx_file.split(' ')[0] + '_' + str(time.time()) + '.tmp'
hgvs_file_name = hgvs_file_name.replace('&','and')
hgvs_file = open(hgvs_file_name, 'w')
successful_local_conversion_results = ''
fail_count = 0
for line in local_conversion_results.split('\n'):
    if 'ERROR' in line or 'None' in line:
        fail_count += 1
        failed_hgvs = line.split('\t')[2]
        hgvs_file.write(failed_hgvs + '\n')
    else:
        successful_local_conversion_results += line

hgvs_file.close()
############
# convert to vcf with VEP
if fail_count > 0:
    vep_call = '/Applications/ensembl-vep/./vep -i ' + hgvs_file_name + ' --port=3337 --refseq --database --force --pick --vcf --output_file STDOUT --no_stats'
    vep_vcf = subprocess.check_output(vep_call, shell = True)
    vep_vcf = vep_vcf.decode('utf-8')
else:
    print("No variants to convert with VEP!")
    vep_vcf = ''
# rm temp file
#subprocess.call(['rm',hgvs_file_name])


# check if anything failed to be included (will happen if HGVS can't be parsed)
# write secondary vcf writing out the missing variants
dat_error_file_name = args.xlsx_file.split(' ')[0] + '_' + panel_type + '.FAILED.dat'
dat_error_file_name = dat_error_file_name.replace('&','and')
dat_errorfile = open(dat_error_file_name, 'w')
for k,v in hgvs_zygosity.items():
	if k not in vep_vcf and k not in successful_local_conversion_results:
		dat_errorfile.write(k + '\t' + str(hgvs_status[k]) + '\n')

# write vcf for gemini annotation, using the zygosity to write the sample genotype
vcf_file_name = args.xlsx_file.split(' ')[0] + '_' + panel_type +  '.vcf'
vcf_file_name = vcf_file_name.replace('&','and')
vcf_file = open(vcf_file_name, 'w')


for line in vep_vcf.split('\n')[:-1]:
    if line[0:2] == '##':
        vcf_file.write(line)
        vcf_file.write('\n')
    elif line[0:6] == '#CHROM':
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=CT,Number=1,Type=String,Description="HGVS to VCF Converter Tool. VEP is VEP, DM is David McGaughey invitae hgvs conversion">\n##FORMAT=<ID=JC,Number=1,Type=String,Description="MVL/John Chiang variant status. Primary is his likely deleterious variants. Secondary is likely benign">\n')
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
	# write variants, with genotype added (0/1, 1/1, etc.)
    else:
        s_line = line.split('\t')
        if line:
			# if conversion fails, then write to error file
            if not s_line[4] or not s_line[5]:
                dat_errorfile.write(s_line[2] + '\t' + str(hgvs_status[s_line[2]]) + '\n')	
            elif 'het' in hgvs_zygosity[s_line[2]]:
                output =  '\t'.join(s_line[0:5]) + '\t' + '.' + '\tPASS\t.\tGT:CT:JC\t' + '0/1:VEP:' + str(hgvs_status[s_line[2]]) + '\n'
                vcf_file.write(output)
            else:
                output =  '\t'.join(s_line[0:5]) + '\t' + '.' + '\tPASS\t.\tGT:CT:JC\t' + '1/1:VEP:' + str(hgvs_status[s_line[2]]) + '\n'
                vcf_file.write(output)

if vep_vcf == '':
	vcf_file.write('##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description="All filters passed">\n##VEP="v89" time="2017-07-28 19:09:25" db="homo_sapiens_core_89_37@ensembldb.ensembl.org" ensembl-funcgen=89.678099a ensembl-io=89.984214b ensembl-variation=89.0afbd49 ensembl=89.ce16adf assembly="GRCh37.p13" gencode="GENCODE 19" genebuild="2011-04" refseq="01_2015"\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=CT,Number=1,Type=String,Description="HGVS to VCF Converter Tool. VEP is VEP, DM is David McGaughey invitae hgvs conversion">\n##FORMAT=<ID=JC,Number=1,Type=String,Description="MVL/John Chiang variant status. Primary is his likely deleterious variants. Secondary is likely benign">\n##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original chr:pos:ref:alt encoding">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
# write the locally converted vcf
for line in local_conversion_results.split('\n')[:-1]:
    s_line = line.split('\t')
    if 'ERROR' not in line and 'None' not in line:
        if 'het' in hgvs_zygosity[s_line[2]]:
            output = line + '\t' + '.' + '\tPASS\t.\tGT:CT:JC\t' + '0/1:DM:' + str(hgvs_status[s_line[2]]) + '\n'
            vcf_file.write(output)
        else:
            output = line + '\t' + '.' + '\tPASS\t.\tGT:CT:JC\t' + '1/1:DM:' + str(hgvs_status[s_line[2]]) + '\n'
            vcf_file.write(output)

dat_errorfile.close()
vcf_file.close()

# vt normalize
subprocess.call(['vt','normalize', vcf_file_name, '-r', '/Users/mcgaugheyd/GenomicData/human_g1k_v37_decoy.fasta', '-o', vcf_file_name[:-4] + '.vt.vcf'])

# sort, bgzip, tabix
subprocess.call(['/Users/mcgaugheyd/git/casey_to_gemini/src/sort_bgzip_tabix.sh', vcf_file_name[:-4] + '.vt.vcf'])

# remove intermediate files
subprocess.call(['rm', vcf_file_name])
subprocess.call(['rm', hgvs_file_name])

print('Done at ' + str(datetime.datetime.now()))
