for line in $(cat $1);
do
 if [ "$(grep -c "^$line" $2)" -ge 1 ]
 then
    grep $line $1 || echo "";
  fi
done
