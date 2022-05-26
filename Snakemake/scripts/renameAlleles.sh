# Edit the output files to match the chromosome numbers in Wragg's data
cut -f1,2,4 OutspeciesInfo_All_aligned.txt | awk '{if ($3 ~ /^\?/) next; else print $0 }' | sed "s/NC_0376//g" | sed "s/\.1/\t/g" | awk -F" " '{print "carnica.LG"$1-37 FS $2 FS $3 FS $4}' > AlignedSnps_focal_CarnicaChr.txt
