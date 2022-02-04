cat EstSfs_Dict0E.csv EstSfs_Dict1E.csv > EstSfs_DictTotal.csv

for cycle in $(seq 2 199)
do
	cat EstSfs_DictTotal.csv EstSfs_Dict${cycle}E.csv > tmp && mv tmp EstSfs_DictTotal.csv
done
