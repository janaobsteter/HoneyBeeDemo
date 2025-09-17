chr=$1
for num in $(seq 0 99)
do
    echo python MapMutations.py ${chr} ${num} 
done
