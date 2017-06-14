#!/usr/local/bin/python3

from openpyxl import load_workbook
import argparse
import subprocess
import datetime

parser = argparse.ArgumentParser(description= 'Takes in xlsx files from John Chiang / MVLGenomics Gene Panels and does the following:\n 1. Identifies hidden rows and skips these \n 2. Extracts hgvs \n\n Usage: casey_xlsx_to_hgvs.py YOURFILE.xlsx > output.tsv')

parser.add_argument('xlsx_file', help = 'John Chiang / OSU / MVL / Casey patient report xlsx file')

args = parser.parse_args()

wb = load_workbook(args.xlsx_file)
ws = wb[wb.get_sheet_names()[0]]
# identify header row, HGVS coding column, and zygosity column
for row in range(2, ws.max_row):
    for column in 'ABCDEFG':
        cell_name = "{}{}".format(column, row)
        if ws[cell_name].value == 'HGVSCoding':
            header_row = row
            HGVS_column = column
        if ws[cell_name].value == 'Zygosity':
            Zygosity_column = column

# get HGVS and zygosity, skipping hidden rows
hgvs_zygosity = []
for row in range(header_row+1, ws.max_row):
    if ws.row_dimensions[row].visible is True:
        cell_name_hgvs = "{}{}".format(HGVS_column, row)
        cell_name_zyg = "{}{}".format(Zygosity_column, row)
        if ws[cell_name_hgvs].value is not None and ws[cell_name_hgvs].value[0:2] == 'NM':
            line = [ws[cell_name_hgvs].value, ws[cell_name_zyg].value]
            hgvs_zygosity.append(line)

# use hgvs_to_vcf.py to convert HGVS to chr, (genomic) coord, ref, alt
for line in hgvs_zygosity:
    print(line)
    vcf = subprocess.check_output(['/Users/mcgaugheyd/git/casey_to_gemini/hgvs_to_vcf.py',line[0]]).decode('utf-8')
    line.append(vcf)
    print(line)


# write failed conversions to error file
# and
# write vcf for gemini annotation
today = datetime.date.today()

error_file_name = args.xlsx_file.split('.xlsx')[0] + '.error'
error_file = open(error_file_name, 'w')
vcf_file_name = args.xlsx_file.split('.xlsx')[0] + '.vcf'
vcf_file = open(vcf_file_name, 'w')
vcf_file.write('##fileformat=VCFv4.2\n')
vcf_file.write('##fileDate=' + str(today).replace('-','') + '\n')
vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')

for line in hgvs_zygosity:
    if line[2].split('\t')[0] == 'ERROR: ':
        error_file.write('\t'.join(line) + '\n')
    else:
        if 'het' in line[1]:
            output = line[2][:-1] + '\t100\tPASS\t.\tGT:GQ:DP\t' + '0/1:100:100\n'
            vcf_file.write(output)
        else:
            output = line[2][:-1] + '\t100\tPASS\t.\tGT:GQ:DP\t' + '1/1:100:100\n'
            vcf_file.write(output)


