# casey_to_gemini
Process to (semi) automate creation of a Gemini database from Casey Eye Institute files

1. On eyeMac (local desktop) run `/git/casey_to_gemini/src/casey_xlsx_to_hgvs.py EXCEL_FILE` which uses the Casey / MVL excel file to parse for HGVS and output a vcf-like file (along with helper files `*FAILED.dat` which lists HGVS variants which couldn't be converted
2. Upload the new vcf and the `*FAILED.dat` files created by `/git/casey_to_gemini/src/casey_xlsx_to_hgvs.py` to `biowulf2:/data/mcgaugheyd/projects/nei/hufnagel/casey_panels/vcf` and `biowulf2:/data/mcgaugheyd/projects/nei/hufnagel/casey_panels/FAILED`, respectively
3. In the `vcf` subfolder run this bash line to create a list of vcfs to process: `for i in *vt.vcf.gz; do wc -l $i; done | grep  -v '^0' | grep -v 'TEMP' |  cut -f2 -d ' ' > all.vcfs`
4. Then create the Gemini database like so: `sbatch --cpus-per-task 16 ~/git/casey_to_gemini/src/create_centralized_gemini_db.py -f all.vcfs 2017-10-16.gemini.db`
5. Write reactive reports with the sbatch a shell script like this:
```
$ cat reports/write5.sh
#!/bin/bash

~/git/casey_to_gemini/write_report_wrapper.sh 2017-11-14.gemini.db 17-2083_RD_v11_assay
~/git/casey_to_gemini/write_report_wrapper.sh 2017-11-14.gemini.db 17-2327_Caroline_
~/git/casey_to_gemini/write_report_wrapper.sh 2017-11-14.gemini.db 17-2440___RD_v11
~/git/casey_to_gemini/write_report_wrapper.sh 2017-11-14.gemini.db 17-2441_v11
```
