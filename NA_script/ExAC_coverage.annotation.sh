
PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
fsource="/home/local/ARCS/nz2274/Resources/ExAC_release/release0.3.1/Coverage"
function annotatedwithMappability {
   fhead="header_add.txt"

   echo "##INFO=<ID=ExAC_covMean,Number=1,Type=String,Description=\"MeancoverageInExAC\">" > $fhead
   echo "##INFO=<ID=ExAC_covMedian,Number=1,Type=String,Description=\"MediancoverageInExAC\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov1,Number=1,Type=String,Description=\"proportion_ExACcoverage>1\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov5,Number=1,Type=String,Description=\"proportion_ExACcoverage>5\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov10,Number=1,Type=String,Description=\"proportion_ExACcoverage>10\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov15,Number=1,Type=String,Description=\"proportion_ExACcoverage>15\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov20,Number=1,Type=String,Description=\"proportion_ExACcoverage>20\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov25,Number=1,Type=String,Description=\"proportion_ExACcoverage>25\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov30,Number=1,Type=String,Description=\"proportion_ExACcoverage>30\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov50,Number=1,Type=String,Description=\"proportion_ExACcoverage>50\">" >> $fhead
   echo "##INFO=<ID=ExAC_cov100,Number=1,Type=String,Description=\"proportion_ExACcoverage>100\">" >> $fhead

   file="/home/local/ARCS/nz2274/Resources/ExAC_release/release0.3.1/Coverage/Panel.ExAC.coverge.bed.gz"
#	for file in $(ls $fsource/Panel.chr*.coverage.txt.gz) 
#		do
   fvcf=$1
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.ExACcov.vcf.gz"
   bcftools annotate -a $file  -c CHROM,FROM,TO,ExAC_covMean,ExAC_covMedian,ExAC_cov1,ExAC_cov5,ExAC_cov10,ExAC_cov15,ExAC_cov20,ExAC_cov25,ExAC_cov30,ExAC_cov50,ExAC_cov100  -h  $fhead  $fvcf |bgzip -c > $fout
#done
#wait
   rm $fhead
}


echo "input format ExAC_coverage.annotation.sh  *.vcf.gz"
#echo " bash $PipePath/ExAC_coverage.annotation.sh  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
annotatedwithMappability   $1 
