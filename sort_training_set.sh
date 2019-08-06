


cd /triton/ics/scratch/users/mpirttin/epigenetics/classification/RFECS/training_set/

for f in K562_enhancers_p300.txt K562_enhancers_p300DNase.txt K562_promoters_DNase.txt K562_random_with_signal.txt K562_promotersDNase_combinedWith_random_with_signal.txt
do
     sort -k1,1 -k2,2n $f > ${f%"."*}"_sorted.txt"

done



