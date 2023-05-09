sort -V AlignedSnps_focal_CarnicaChr.txt | awk -F" " '{print $0 FS $1"_"$2}' > AlignedSnps_focal_CarnicaChrSorted.txt
grep -Fwf VcfFullPos.txt AlignedSnps_focal_CarnicaChrSorted.txt > AlignedSnps_focal_CarnicaChrSortedVcf.txt
