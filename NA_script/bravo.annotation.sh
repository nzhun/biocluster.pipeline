
PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
function annotatedBRAVO {
  
   fMAP=$1
   fvcf=$2
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.BRAVO_MAF.vcf.gz"
   echo $fMAP

#  if [ ! -f $fMAP ]; then
 #        Redo_mappability $len
 #  fi
   keywd="BRAVO_AF"
   des="Bravo population frequency in file $fMAP"
   fhead="header_add.txt"

   echo "##INFO=<ID=$keywd,Number=1,Type=String,Description=\"$des\">" > $fhead
   bcftools annotate -a $fMAP  -c CHROM,FROM,TO,REF,ALT,BRAVO_AF  -h  $fhead  $fvcf |bgzip -c > $fout
}


echo "input format annotatedBRAVO  *bed.gz  *.vcf.gz"
echo " example: bash $PipePath/annotatedBRAVO.sh /home/local/ARCS/nz2274/Resources/BRAVO/bravo.topmed.bed.gz  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
echo $#
if [[ $# -lt 2 ]]; then
    exit;
fi

#source $PipePath/mappability_generate.sh

annotatedBRAVO   $1 $2 
if [ ! -f $1 ];then
   echo "$1 cannot find!\n"
   exit;
fi

if [ ! -f $2 ]; then
  echo "$2 cannot find!\n"
  exit;
fi
