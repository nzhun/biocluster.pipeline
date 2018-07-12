
PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
fsource="/home/local/users/jw/resources/gnomad/"
function annotatedwithcoverage {
   fvcf=$1
   for type in "Genome" # "Exome"
   do	
	   fhead="header_add.txt"
	   echo "##INFO=<ID=GnomAD_${type}_covMean,Number=1,Type=String,Description=\"MeancoverageInGnom_${type}_AD\">" > $fhead
	   echo "##INFO=<ID=GnomAD_${type}_covMedian,Number=1,Type=String,Description=\"MediancoverageInGnom_${type}_AD\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov1,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>1\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov5,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>5\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov10,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>10\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov15,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>15\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov20,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>20\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov25,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>25\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov30,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>30\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov50,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>50\">" >> $fhead
	   echo "##INFO=<ID=GnomAD_${type}_cov100,Number=1,Type=String,Description=\"proportion_GnomAD_${type}_coverage>100\">" >> $fhead

  
	 
			file="$fsource/$type/$type.cov.txt.gz"		
				fbname=$(basename $fvcf)
				fbname=${fbname%.vcf.gz}
				fbname=${fbname%.vcf}
				fout="${fbname}.GnomADcov.vcf.gz"
				bcftools annotate -a $file  -c CHROM,POS,GnomAD_${type}_covMean,GnomAD_${type}_covMedian,GnomAD_${type}_cov1,GnomAD_${type}_cov5,GnomAD_${type}_cov10,GnomAD_${type}_cov15,GnomAD_${type}_cov20,GnomAD_${type}_cov25,GnomAD_${type}_cov30,GnomAD_${type}_cov50,GnomAD_${type}_cov100  -h  $fhead  $fvcf |bgzip -c > $fout
		
	done

}


echo "input format GnomAD_coverage.annotation.sh  *.vcf.gz"
#echo " bash $PipePath/GnomAD_coverage.annotation.sh  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
annotatedwithcoverage   $1
if [ ! -f $1 ] ; then
   echo $1" cannot find!"
fi 
