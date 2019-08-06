


cd /triton/ics/project/csb/software/public/RFECS/K562_data

for f in `ls`; do
sort -k1,1 -k2,2n $f > ${f%"."*}"sorted.bed"

done


cd /triton/ics/project/csb/software/public/RFECS/GM12878_data

for f in `ls`; do
sort -k1,1 -k2,2n $f > ${f%"."*}"sorted.bed"

done
