 
#for i in {1..22} X Y
#do
#	nm="gnomad.genomes.chr$i.knowngenes"
	#tabix $gnomeAD/Genome/gnomad.genomes.r2.0.1.sites.$i.vcf.gz -R ~/PAH/PAH-CHD/source/CHD_PAH.knownGenes.bed -h |bgzip -c  > $nm.vcf.gz
 #   perl ~/Pipeline/NA_script/vcf2bed_gnomAD.pl $nm.vcf.gz &
	#bedtools intersect -a $nm.vcf.gz.2.bed -b $XGEN -wa -header > $nm.vcf.gz.xgen.bed
#	file="$nm.vcf.gz.2.bed"
	#awk 'BEGIN{FS="\t";OFS="\t"}{if($1 ~/^#/){print;next;}a=$4;b=$5;if(length($4)>1 && length($4)<length($5)) {a=substr($4,0,1);b=substr($5,0,length($5)-length($4)+1)}if(length($5)>1 && length($5)<=length($4)) {a=substr($5,0,1);b=substr($4,0,length($4)-length($5)+1);};$4=a;$5=b;$2=$2+1;$3=$2+length($4)-1;print }' $bedf > $nm.corrected.bed
	#file="$nm.corrected.bed"
	#perl $ANNOVAR/table_annovar.pl $file $ANNHDB --buildver hg19 --remove -protocol revel,mcap,dbnsfp33a,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,exac03,gnomad_exome,gnomad_genome -operation f,f,f,f,f,f,f,f,f -otherinfo  -nastring . 
#done


##get the NFE number 

## restrict to NFE PAH-CHD
file=$1

perl $ANNOVAR/table_annovar.pl $file $ANNHDB --buildver hg19 --remove -protocol refGene -operation g -otherinfo  -nastring . $2
