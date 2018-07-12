
PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
function annotatedwithMappability {
   len=$1
   fMAP=$2
   fvcf=$3
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.mappability.vcf.gz"
   echo $fMAP
   fbm=$(basename $fMAP)
   pref=$(echo $fbm|cut -f 1 -d"_")
   
  if [ ! -f $fMAP ]; then
         Redo_mappability $len	$REF $pref
   fi
   keywd="Mappability"
   des="Mappability is obtained in file $fMAP"
   fhead="header_add.txt"

   echo "##INFO=<ID=$keywd,Number=1,Type=String,Description=\"$des\">" > $fhead
   bcftools annotate -a $fMAP  -c CHROM,FROM,TO,Mappability  -h  $fhead  $fvcf |bgzip -c > $fout
}


echo "input format annotatedwithMappability len *bed.gz  *.vcf.gz"
echo " example: bash $PipePath/mappability.annotation.sh 152 /home/local/ARCS/nz2274/Resources/mappability/hg19_152bp_mappability.bed.gz  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
echo $#
if [[ $# -lt 3 ]]; then
    exit;
fi

source $PipePath/mappability_generate.sh 

annotatedwithMappability   $1 $2 $3
if [ ! -f $2 ];then
   echo "$2 cannot find!\n"
fi

if [ ! -f $3 ]; then
  echo "$3 cannot find!\n"
fi
