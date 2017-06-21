#!/bin/bash

vcf=$1

grep '^#' "$vcf" > "$vcf".temp
grep -v '^#' "$vcf" | sort -k1,1 -k2,2n >> "$vcf".temp
mv "$vcf".temp "$vcf"
bgzip "$vcf"
tabix -p vcf "$vcf".gz
