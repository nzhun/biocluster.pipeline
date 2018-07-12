


## download  picard.jar   http://sourceforge.net/projects/picard/


  cat ~/solexa/home/petkraw/References/human_g1k_v37.dict> clinvar.interval.list
  awk 'BEGIN{FS="\t"}{OFS="\t"; if($1=="-"||$2=="-"||$3=="-"){next;} if($3!="-"&&$4!="-"){next} print  $1,$2-(length($3)*15),($2+length($3)-1),"+\t"$1","($2-(length($3)*15))","$10","$3","$4; }' ../../projects/eclipse_14.04/XMLParser/clinvar.csv |sort -k1,1n -k2,2n >> clinvar.interval.list

 java -jar picard.jar  ExtractSequences  INTERVAL_LIST=clinvar.interval.list  R=~/solexa/home/petkraw/References/human_g1k_v37.fasta  O=out.list2


cp out.list2 out.list
 awk 'BEGIN{ls="";}{if($_ ~/-$/){if(ls !=""){print ls;}ls=$_",";next}else{ls=ls""$_}}' out.list >temp
mv temp out.list

tr -d '\n' <out.list|sed -e 's/>/\n>/g' |grep -v '^$'|sed -e 's/>//g'|awk 'BEGIN{OFS="\t";FS=",";}{ l=length($4);al=length($6);for(j=al;j>0;j--){ if(substr($6,j-l,j)==$4){}} print $1,$2,".",$6,substr($6,0, length($6)-length($4)+1),"100","PASS","AC="$3,"GT","0/1";}' >clinvar.indel.vcf


#java -jar ../../../Applications/GATK3.3/GenomeAnalysisTK.jar   -R ~/solexa/home/petkraw/References/human_g1k_v37.fasta  -T LeftAlignAndTrimVariants -trim
# --variant clinvar.indel.vcf  -o norm.vcf

#trim file
awk '{OFS="\t";if($1 ~/^#/){print $_;next} al=length($4);l=length($4)-length($5); str=substr($4,al-l+1,l);st=al-l; astr=""; for(j=al;j>0;j=j-l){ st1=substr($4,j-l+1,l);astr=astr","st1; if(st1!=str){st=j;break;}} print $1,$2+st-1,$3,substr($4,st,l+1),substr($4,st,1),$6,$7,$8 }' clinvar.indel.vcf|less


awk '{OFS="\t";if($1 ~/^#/){print $_;next} al=length($4);l=length($4)-length($5); str=substr($4,al-l+1,l);st=al-l; astr=""; for(j=al;j>0;j=j-l){ st1=substr($4,j-l+1,l);astr=astr","st1; if(st1!=str){st=j;break;}} print $1,$2+st-1,$3,substr($4,st,l+1),substr($4,st,1),$6,$7,$8 }' clinvar.indel.vcf >correct.indel.corrected.vcf


perl ../../projects/eclipse_14.04/XMLParser/upinfo.pl clinvar.indel.corrected.vcf  clinvar.csv clinvar.corrected.csv



