vcf=$1;
nm=$(basename $vcf)
mkdir TEMP
nm=${nm%.vcf*}
out_plink=$nm".plink" 
if [[ "$vcf" == *.vcf.gz ]]; then
{
#	nm=${nm%.vcf.gz}
#	out_plink=$nm".plink"
	if [ ! -f $out_plink".ped" ]; then	
                cmd="vcftools --gzvcf $vcf  --out $out_plink --temp  TEMP --maf 0.05  --plink"
                echo $cmd
		        vcftools --gzvcf $vcf --plink-tped --out $out_plink.tmp --temp  TEMP --maf 0.05 
		        plink --tfile $out_plink.tmp --out $out_plink
		#vcftools --gzvcf $vcf  --out $out_plink --temp  TEMP --maf 0.05  --plink  # --relatedness  #--relatedness2  #
	fi
}
else
{
	
	if [ ! -f $out_plink".ped" ]; then	
		  vcftools --gzvcf $vcf --plink-tped --out $out_plink.tmp --temp  TEMP --maf 0.05 # vcftools --vcf $vcf --out $out_plink --temp TEMP --maf 0.05 --plink  # --relatedness  #--relatedness2  #
	fi
}	
	
fi
#outplink_hapmap="HAPMAP.plink"
hapmap="/home/local/ARCS/nz2274/Resources/hapmap_r23a/hapmap_r23a"
#vcftools --vcf $HAPMAP --out $outplink_hapmap --temp TEMP --maf 0.05 --plink  # --relatedness  #--relatedness2  #
mergedplink=$nm"_merged.plink"

plink --bfile $out_plink --bmerge $hapmap.bed $hapmap.bim $hapmap.fam --recode --out $mergedplink

#excludes=$(awk 'BEGIN{str=""}{str=str","$1}END{str=substr(str,1,length(str)-1); print str}' $mergedplink.missnp)
 
ret=$(plink --file $out_plink --exclude  $mergedplink.missnp --out $out_plink"_correct" --make-bed)
plink --bfile $out_plink"_correct" --bmerge $hapmap.bed $hapmap.bim $hapmap.fam --recode --out $mergedplink

echo "retrun is ".$ret
echo $out_plink

out=$nm"_hapmap.relationships"

plink --file $mergedplink --genome --min 0.05 --out $out 
