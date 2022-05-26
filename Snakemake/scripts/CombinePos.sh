chr=( "$@" )
files=()
for val in ${chr[@]}
do
        files+=(FocalChrPos$val.txt)
done

cat "${files[@]}" > FocalChrPos_All.txt
# Remove the lines without the sites
grep -v Amel FocalChrPos_All.txt > tmp
grep -v Acer tmp > tmp1
grep -v Ador tmp1 > tmp2
grep -v Aflo tmp2 > tmp3
mv tmp3 FocalChrPos_All.txt
rm tmp tmp1 tmp2
