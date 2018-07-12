
#PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
fsource="/home/local/users/jw/resources/gnomad/"

function generate_header {
  less  $fsource/Genome/gnomad.genomes.r2.0.1.sites.1.vcf.gz|head -n 300|grep '##INFO'|egrep 'AC,|AF,|AN,|AC_AFR,|AC_AMR,|AC_ASJ,|AC_EAS,|AC_FIN,|AC_NFE,|AC_OTH,|AC_SAS,|AC_Male,|AC_Female,|AN_AFR,|AN_AMR,|AN_ASJ,|AN_EAS,|AN_FIN,|AN_NFE,|AN_OTH,|AN_SAS,|AN_Male,|AN_Female,|AF_AFR,|AF_AMR,|AF_ASJ,|AF_EAS,|AF_FIN,|AF_NFE,|AF_OTH,|AF_SAS,|AF_Male,|AF_Female,|POPMAX,|AC_POPMAX,|AN_POPMAX,|AF_POPMAX,|DP_MEDIAN,|DREF_MEDIAN'  |sed 's/ID=/ID=gnomAD_Genome_/g' > /home/local/users/jw/resources/gnomad/Genome/header.txt
  
  less  $fsource/Exome/gnomad.exomes.r2.0.1.sites.vcf.gz|head -n 300|grep '##INFO'|egrep 'AC,|AF,|AN,|AC_AFR,|AC_AMR,|AC_ASJ,|AC_EAS,|AC_FIN,|AC_NFE,|AC_OTH,|AC_SAS,|AC_Male,|AC_Female,|AN_AFR,|AN_AMR,|AN_ASJ,|AN_EAS,|AN_FIN,|AN_NFE,|AN_OTH,|AN_SAS,|AN_Male,|AN_Female,|AF_AFR,|AF_AMR,|AF_ASJ,|AF_EAS,|AF_FIN,|AF_NFE,|AF_OTH,|AF_SAS,|AF_Male,|AF_Female,|POPMAX,|AC_POPMAX,|AN_POPMAX,|AF_POPMAX,|DP_MEDIAN,|DREF_MEDIAN'  |sed 's/ID=/ID=gnomAD_Exome_/g' > /home/local/users/jw/resources/gnomad/Exome/header.txt
  

}
function getBaseName {
	file=$1
	fbname=$(basename $fvcf)
	fbname=${fbname%.vcf.gz}
	fbname=${fbname%.vcf}
	echo $fbname
}


function callbcftools {
	file=$1 
	type=$2
	fhead=$3
	fvcf=$4
	fout=$5	 	cols="CHROM,FROM,TO,REF,ALT,gnomAD_${type}_AC,gnomAD_${type}_AF,gnomAD_${type}_AN,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,gnomAD_${type}_AC_AFR,gnomAD_${type}_AC_AMR,gnomAD_${type}_AC_ASJ,gnomAD_${type}_AC_EAS,gnomAD_${type}_AC_FIN,gnomAD_${type}_AC_NFE,gnomAD_${type}_AC_OTH,gnomAD_${type}_AC_SAS,gnomAD_${type}_AC_Male,gnomAD_${type}_AC_Female,gnomAD_${type}_AN_AFR,gnomAD_${type}_AN_AMR,gnomAD_${type}_AN_ASJ,gnomAD_${type}_AN_EAS,gnomAD_${type}_AN_FIN,gnomAD_${type}_AN_NFE,gnomAD_${type}_AN_OTH,gnomAD_${type}_AN_SAS,gnomAD_${type}_AN_Male,gnomAD_${type}_AN_Female,gnomAD_${type}_AF_AFR,gnomAD_${type}_AF_AMR,gnomAD_${type}_AF_ASJ,gnomAD_${type}_AF_EAS,gnomAD_${type}_AF_FIN,gnomAD_${type}_AF_NFE,gnomAD_${type}_AF_OTH,gnomAD_${type}_AF_SAS,gnomAD_${type}_AF_Male,gnomAD_${type}_AF_Female,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,gnomAD_${type}_POPMAX,gnomAD_${type}_AC_POPMAX,gnomAD_${type}_AN_POPMAX,gnomAD_${type}_AF_POPMAX,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,- "
	if [ "$type" == "Genome" ]; then
	cols="CHROM,FROM,TO,REF,ALT,gnomAD_Genome_AC,gnomAD_Genome_AF,gnomAD_Genome_AN,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,gnomAD_Genome_AC_AFR,gnomAD_Genome_AC_AMR,gnomAD_Genome_AC_ASJ,gnomAD_Genome_AC_EAS,gnomAD_Genome_AC_FIN,gnomAD_Genome_AC_NFE,gnomAD_Genome_AC_OTH,gnomAD_Genome_AC_Male,gnomAD_Genome_AC_Female,gnomAD_Genome_AN_AFR,gnomAD_Genome_AN_AMR,gnomAD_Genome_AN_ASJ,gnomAD_Genome_AN_EAS,gnomAD_Genome_AN_FIN,gnomAD_Genome_AN_NFE,gnomAD_Genome_AN_OTH,gnomAD_Genome_AN_Male,gnomAD_Genome_AN_Female,gnomAD_Genome_AF_AFR,gnomAD_Genome_AF_AMR,gnomAD_Genome_AF_ASJ,gnomAD_Genome_AF_EAS,gnomAD_Genome_AF_FIN,gnomAD_Genome_AF_NFE,gnomAD_Genome_AF_OTH,gnomAD_Genome_AF_Male,gnomAD_Genome_AF_Female,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,gnomAD_Genome_POPMAX,gnomAD_Genome_AC_POPMAX,gnomAD_Genome_AN_POPMAX,gnomAD_Genome_AF_POPMAX,-,-,-,-,-,-,-,-,-"
	fi
	echo "bcftools annotate -a $file  -c $cols -h  $fhead  $fvcf | bgzip -c  > $fout"	
	bcftools annotate -a $file  -c $cols -h  $fhead  $fvcf | bgzip -c  > $fout 
	#echo $cmd 
	#$cmd
}

function separatevcf {
	file=$1
	baseorg=$(getBaseName $file)
	for chr in {1..22} X Y
	do
		tabix $file $chr -h |bgzip -c > $baseorg.$chr.vcf.gz
		tabix -f -p vcf $baseorg.$chr.vcf.gz
	done
    echo $baseorg
}
function mergevcf {
	str=$1;
}
function annotatedwithPopulation {
    fvcf=$1
	type="Exome"
    fhead="$fsource/$type/header.txt"  #header_add.txt
	if [ ! -f $fhead ]; then
		generate_header
	fi
	
	file="$fsource/$type/gnomad.exomes.r2.0.1.sites.vcf.gz.2.bed.gz"		
	fbname=$(getBaseName $file)
	fout="${fbname}.GnomAD_${type}_POP.vcf.gz"
	callbcftools $file $type $fhead $fvcf $fout
	tabix -f -p vcf $fout
	
	## genome
	
	fprefix=$(separatevcf $fout)
	type="Genome"
    fhead="$fsource/$type/header.txt"  #header_add.txt
	str_merge="vcf-concat ";
	for chr in {1..22} X Y 
 		do
			{
			fvcf_sub=$fprefix.$chr.vcf.gz
			file="$fsource/$type/gnomad.genomes.r2.0.1.sites.$chr.vcf.gz.2.bed.gz"		
 			fbname=$(getBaseName $fvcf_sub)
 			fout="${fbname}.GnomAD_${type}_POP.$chr.vcf.gz"
			echo "callbcftools $file $type $fhead $fvcf_sub $fout"
 			callbcftools $file $type $fhead $fvcf_sub $fout
			str_merge=$str_merge" "$fout
		} &
		done
	wait
 	
	#str_merge=$str_merge" |bgzip -c > ${fbname}.GnomAD_${type}_POP.vcf.gz "
	#echo $str_merge
	#$str_merge
}


echo "input format gnomeAD_population.sh  *.vcf.gz"
#echo " bash $PipePath/GnomAD_coverage.annotation.sh  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
if [  -z "$1" ];then
  echo "vcf input is required!";
  exit;
fi
annotatedwithPopulation   $1 
