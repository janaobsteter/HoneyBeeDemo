#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

HEADERNUM=$(( $(grep "##" $1 | wc -l) + 1 ))
INFOLINE=$(( $(grep -Fn "INFO" $2 | cut --delimiter=":" --fields=1  | head -n1) ))
awk -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<HEADER{{{{print}}}}; \
FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"}}}}; \
FNR>HEADER{{{{$8=a[FNR-HEADER]; print}}}}' OFS="\t" $1 $2
