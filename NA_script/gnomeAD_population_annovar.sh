
#PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
fsource="/home/local/users/jw/software_packages/annovar/humandb/"

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
	fout=$5	 	
	cols="CHROM,FROM,TO,REF,ALT,gnomAD_Exome_ALL,gnomAD_Exome_AF_AFR,gnomAD_Exome_AF_AMR,gnomAD_Exome_AF_ASJ,gnomAD_Exome_AF_EAS,gnomAD_Exome_AF_FIN,gnomAD_Exome_AF_NFE,gnomAD_Exome_AF_OTH,gnomAD_Exome_AF_SAS"
	if [ "$type" == "Genome" ]; then
	cols="CHROM,FROM,TO,REF,ALT,gnomAD_Genome_ALL,gnomAD_Genome_AF_AFR,gnomAD_Genome_AF_AMR,gnomAD_Genome_AF_ASJ,gnomAD_Genome_AF_EAS,gnomAD_Genome_AF_FIN,gnomAD_Genome_AF_NFE,gnomAD_Genome_AF_OTH"
	fi
		
	bcftools annotate -a $file  -c $cols -h  $fhead  $fvcf | bgzip -c  > $fout 
	#echo $cmd 
#	$cmd
}


function annotatedwithPopulation {
    fvcf=$1
	type="Exome"
    fhead="$fsource/gnomAD.$type.header.txt"  #header_add.txt
	if [ ! -f $fhead ]; then
		generate_header
	fi
	
	file="$fsource/hg19_gnomad_exome.bed.gz"		
	fbname=$(getBaseName $file)
	fout="${fbname}.GnomAD_${type}_POP.vcf.gz"
	callbcftools $file $type $fhead $fvcf $fout
	tabix -f -p vcf $fout
	
	## genome
	
	#fprefix=$(separatevcf $fout)
	type="Genome"
    fhead="$fsource/gnomAD.$type.header.txt"  #header_add.txt
	
			
	file="$fsource/hg19_gnomad_genome.bed.gz"		
 	fbname=$(getBaseName $fvcf)
 	fout_genome="${fbname}.GnomAD_${type}_POP.vcf.gz"
	echo "callbcftools $file $type $fhead $fvcf_sub $fout"
 	callbcftools $file $type $fhead $fout $fout_genome
	
 	
}


echo "input format gnomeAD_population.sh  *.vcf.gz"
#echo " bash $PipePath/GnomAD_coverage.annotation.sh  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
annotatedwithPopulation   $1 
