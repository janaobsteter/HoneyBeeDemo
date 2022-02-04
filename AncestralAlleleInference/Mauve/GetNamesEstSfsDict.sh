for cycle in $(seq 0 199)
do
	cut -f1 EstSfs_Dict${cycle}.csv > EstSfsNames${cycle}.csv
done
