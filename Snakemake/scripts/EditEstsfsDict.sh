for cycle in $(seq 999)
do
  	cut -f2,3,4,5 EstSfs_Dict${cycle}.csv > EstSfs_Dict${cycle}E.csv
        grep -v "()" EstSfs_Dict${cycle}E.csv > tmp && mv tmp EstSfs_Dict${cycle}E.csv
        sed -i "s/ //g" EstSfs_Dict${cycle}E.csv
        # Set the correct separators for the file
        awk -F "\t" '{print $1"\t"$2" "$3" "$4}' EstSfs_Dict${cycle}E.csv > tmp && mv tmp EstSfs_Dict${cycle}E.csv
        # Remove the parenthesis from the file
        sed -i "s/(//g" EstSfs_Dict${cycle}E.csv
        sed -i "s/)/ /g" EstSfs_Dict${cycle}E.csv
done
