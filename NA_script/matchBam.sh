## match vcf id to bam file
fbam="../../../BAM/bam.txt"  ### format: bam file path per line
fout="variants.bam.txt" ###
fvariants="variants.txt" ### format, chr\tpot\tind1|ind2...
#cline=$(wc -l variants.txt |cut -f 1 -d" ")
rm $fout
while read -r line
do
	query=$(echo $line|cut -f 3)
	rs=$(egrep "$query" $fbam|awk 'BEGIN{str=""}{m=split($1,ph,"/");$1=ph[m]; if(str==""){str=$1}else{str=str","$1}}END{print str}')
	echo -e $line"\t"$rs|sed -e 's/ /	/g'|cut -f 1-2,4 >>$fout
done < "$fvariants"