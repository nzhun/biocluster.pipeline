out=$1
head -n 1 $out|sed  's/^/#/g' > $out.temp.bed
egrep '^[1-9]' $out|sort -k1,1n -k2,2n -k3,3n >>$out.temp.bed
egrep '^[X]' $out|sort -k1,1n -k2,2n -k3,3n  >>$out.temp.bed
egrep '^[Y]' $out|sort -k1,1n -k2,2n -k3,3n >>$out.temp.bed

echo "output is in $out.temp.bed"
