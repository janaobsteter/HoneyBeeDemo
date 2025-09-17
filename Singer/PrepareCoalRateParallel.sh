for chromosome in 14 15 16
do
    echo python GetCoalescenceRates_Lineages.py ${chromosome}
    echo python GetCoalescenceRates_Subspecies.py ${chromosome}
    echo python GetCoalescenceRates_Country.py ${chromosome}
done

