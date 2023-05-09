chr=( "$@" )
files=()
for val in ${chr[@]}
do
        files+=(OutspeciesInfo$val.txt)
done

cat "${files[@]}" > OutspeciesInfo_All.txt
grep -v -w "?,?,?" OutspeciesInfo_All.txt > OutspeciesInfo_All_aligned.txt
awk -F"\t" '{print $1 FS $2}' OutspeciesInfo_All_aligned.txt > FocalChrPos_All_aligned.txt
