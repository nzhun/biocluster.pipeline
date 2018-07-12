file=$1  ### vcf tabixed 
bed=$2; ## bed format
n=$3  ## threads to use
interval=$4 ## length for single thread
echo "totathreads : "$n
echo "interval: $interval"
### target of this script is to split $file into n parts and processing it parrallel
lc=$(wc $bed|awk '{print $1}')
i=0;
for (( j=$interval ; j <= $lc + $interval ; j = j + $interval ))
do
   i=$((i+1))
   if [ -f "$file.temp.$i.vcf" ]; then
	   continue;
   fi
   head -n $j $bed|tail -n $interval >temp.$i.bed
   opt="-h"
 
 
   echo "tabix $file -R temp.$i.bed $opt  > $file.temp.$i.vcf"
   if [ $(( i % $n )) == 0 ]; then
{
 
  	tabix $file -R temp.$i.bed $opt   > $file.temp.$i.vcf
    cat $file.temp.$i.vcf|vcf-sort -c |bgzip -c >$file.$i.vcf.gz
	if [ $? == 0 ]; then 
	   tabix $file.$i.vcf.gz
	   rm  $file.temp.$i.vcf
	fi
   
}
   else 
    {
	  	tabix $file -R temp.$i.bed $opt   > $file.temp.$i.vcf
	    cat $file.temp.$i.vcf|vcf-sort -c |bgzip -c >$file.$i.vcf.gz
		if [ $? == 0 ]; then 
		   tabix $file.$i.vcf.gz
		   rm  $file.temp.$i.vcf
		fi
   }&
	fi
   
done


wait
echo "output $i parts "


dst="$file.subset.merged.vcf"
str=""
rm $dst
for (( k=1; k<= $i;k++ ))
do
	str=$str" "$file.temp.$i.vcf.gz
done
bgzip $dst
tabix $dst.gz
	
	

