#########################################################################
# File Name: ../00_Codes/intersect_split.sh
# Author: chuan jiang 
# mail: chujiang@ttuhsc.edu
# Created Time: Tue Sep 23 10:36:01 2025
#########################################################################
#!/bin/bash
# split_bed.sh
# Usage: bedtools intersect ... | ./split_bed.sh overlap_file nonoverlap_file

# 参数
OVERLAP_FILE="$1"
NONOVERLAP_FILE="$2"

if [[ -z "$OVERLAP_FILE" || -z "$NONOVERLAP_FILE" ]]; then
    echo "Usage: $0 <overlap_file> <nonoverlap_file>"
    exit 1
fi

awk -v overlap="$OVERLAP_FILE" -v nonoverlap="$NONOVERLAP_FILE" '
{
    if ($NF == "." || $NF == 0) 
        print > nonoverlap
    else 
        print > overlap
}
' 

