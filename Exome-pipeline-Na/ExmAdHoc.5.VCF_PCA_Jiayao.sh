#!/bin/bash
###############################################################
# PCA for population structure
# Jiayao Wang 0815 2017
###############################################################
usage="
ExmAdHoc.5.VCF_PCA.sh -i <InputFile> -l <logfile> -H
     -i (required) - A vcf input file
     -o (optional) - output name - will be derived from input filename if not provided
     -H (flag) - echo this message and exit
"

SamOnly="false"
while getopts i:r:o:b:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        o) OutNam="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#some variables
EXOMFILT=/home/local/users/jw/CUMC/Exome-Filters-Jiayao/Filter_For_PCA.py
source $RefFil
echo $HapMapReference


InpFil=`readlink -f $InpFil`
if [[ -z "$OutNam" ]];then OutNam=`basename $InpFil`; OutNam=${OutNam/.bed/}; OutNam=${OutNam/.vcf/}; fi # a name for the output files

LogFil=$OutNam.PCA.log

#check for vcf, if vcf convert to plink format
if [[ "${InpFil##*.}" != "bed" ]]; then
    VCFFil=`readlink -f $InpFil`
    BbfNam=$OutNam
    #filter the VCF for common variants by rsid
    #$EXOMFILT/ExmFilt.0.FilterbyrsID.py -v $VcfFil -o $BbfNam
   if [ ! -f $BbfNam.filter.aaf.vcf.gz  ]; then
	   FilteredOut=$(basename ${VCFFil}).filtered.vcf
	   echo `pwd`
	   echo $EXOMFILT -v $VCFFil -o $FilteredOut -c $HapMapReference
	   $EXOMFILT -v $VCFFil -o $FilteredOut -c $HapMapReference
	   mkdir sort_$(basename ${VCFFil})
	   vcf-sort $FilteredOut  -t ./sort_$(basename ${VCFFil}) > $BbfNam.filter.aaf.vcf
    VcfFil=$BbfNam.filter.aaf.vcf
    bgzip $BbfNam.filter.aaf.vcf
    VcfFil=$BbfNam.filter.aaf.vcf.gz
    tabix -p vcf $VcfFil
	rm -rf $FilteredOut sort_$(basename ${VCFFil})
	else
		echo "Filtered VCF present"
		VcfFil=$BbfNam.filter.aaf.vcf.gz	
    fi    
# convert vcf --> plink using vcftools; if there are more than 1000 samples in the vcf vcftools cannot generate the necesary number of temporary files (cluster limitiation) so we need to do it in multiple batches and then remerge
	touch ExcludeSNPs.list
    SamLst=$OutNam.samples.list
    less $VcfFil | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' >> $SamLst
    SamLen=`cat $SamLst | wc -l`
    if [[ $SamLen -lt 500 ]]; then
		echo vcftools --gzvcf $VcfFil --plink --out $BbfNam.tmp
        vcftools --gzvcf $VcfFil --plink --out $BbfNam.tmp
		echo plink --file $BbfNam.tmp --make-bed --out $BbfNam
        plink --file $BbfNam.tmp --make-bed --out $BbfNam
        rm -f $BbfNam.tmp*
    else
		echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		echo "PLINK THE CASE"
		echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        mkdir -p TempSplit
        split --additional-suffix=.list -l 500 $SamLst $OutNam.samples.split.
        for i in $OutNam.samples.split.*.list; do
			echo $i
            echo vcftools --gzvcf $VcfFil --keep $i --exclude ExcludeSNPs.list --plink --out TempSplit/${i/.list/} &
            vcftools --gzvcf $VcfFil --keep $i --exclude ExcludeSNPs.list --plink --out TempSplit/${i/.list/} &
            #vcftools --gzvcf $VcfFil --keep $i --plink --out TempSplit/${i/.list/}
            echo plink --file TempSplit/${i/.list/} --make-bed --out TempSplit/${i/.list/}
            plink --file TempSplit/${i/.list/} --make-bed --out TempSplit/${i/.list/}
        done
        wait
		echo $BbfNam.tmp.splitlist
        ls TempSplit/*ped | awk -v OFS='\t' '{ $2=$1 ; gsub ( /ped$/, "map", $2) ; print }' > $BbfNam.tmp.splitlist
		echo plink --merge-list $BbfNam.tmp.splitlist --make-bed --out $BbfNam
        plink --merge-list $BbfNam.tmp.splitlist --make-bed --out $BbfNam
        rm -rf TempSplit
    fi
    
    echo
    echo "------------------------------------------------------------------------"
    echo
    #change -9 in the fam to 2
    awk '{ gsub( /-9$/, "2");if(length($1)>15){$1=substr($1,length($1)-15);$2=substr($2,length($2)-15)} print }' $BbfNam.fam > $BbfNam.fam.temp
    mv -f $BbfNam.fam.temp $BbfNam.fam
else
    BbfNam=`readlink -f $InpFil`
    BbfNam=${BbfNam/.bed/}
    echo $BbfNam
    awk '{ gsub( /-9$/, "2"); print }' $BbfNam.fam > $BbfNam.fam.temp
    mv $BbfNam.fam $BbfNam.fam.pcabkp
    mv $BbfNam.fam.temp $BbfNam.fam
fi
folder=$(dirname $BbfNam)
# remove LD
plink --bfile $BbfNam --indep-pairwise 50 5 0.5
SnpList=plink.prune.in


# Get HapMap data
bash /home/local/users/jw/resources/AncestryPCA/PLINK/MakePlink.sh 1KG.vcf
HapMapDat=1KG.XGEN.SNP.Common

#HapMapDat=$HapMapReference #$OutNam"_HapMapData"
echo $HapMapDat
if [[ $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo


EigDat=$OutNam.HapMap
plink --bfile $BbfNam --bmerge $HapMapDat.bed $HapMapDat.bim $HapMapDat.fam --geno 0.0005 --allow-no-sex --make-bed --out $EigDat

#check for mismatched snps and multiple position/chr snps and exclude and remerge if necessary
grep "Warning: Multiple [cp]" $EigDat.log | sed s/.*\ \'//g | sed s/\'.*//g > ExcludeSNPs.list
cat $EigDat-merge.missnp >> ExcludeSNPs.list

if [[ -s ExcludeSNPs.list ]]; then
    plink --bfile $HapMapDat --exclude ExcludeSNPs.list --allow-no-sex --make-bed --out $HapMapDat
    plink --bfile $BbfNam --exclude ExcludeSNPs.list --allow-no-sex --make-bed --out $BbfNam  

    if [[ $? -ne 0 ]]; then exit; fi
    echo
    echo "------------------------------------------------------------------------"
    echo
    plink --bfile $BbfNam --bmerge $HapMapDat.bed $HapMapDat.bim $HapMapDat.fam --geno 0.0005 --allow-no-sex --make-bed --out $EigDat
    if [[ ! -e $EigDat-merge.missnp && $? -ne 0 ]]; then exit; fi
    echo
    echo "------------------------------------------------------------------------"
    echo
fi
awk '{ gsub( /-9$/, "1"); print }' $EigDat.fam > $EigDat.fam.temp
mv -f $EigDat.fam.temp $EigDat.fam


# Convert data to Eigenstrat format
cp $EigDat.fam $EigDat.pedind
echo genotypename: $EigDat.bed > par.BBF.EIGENSTRAT
echo snpname: $EigDat.bim >> par.BBF.EIGENSTRAT
echo indivname: $EigDat.pedind >> par.BBF.EIGENSTRAT
echo outputformat: EIGENSTRAT >> par.BBF.EIGENSTRAT
echo genooutfilename: $OutNam.eigenstratgeno >> par.BBF.EIGENSTRAT
echo snpoutfilename: $OutNam.snp >> par.BBF.EIGENSTRAT
echo indoutfilename: $OutNam.ind >> par.BBF.EIGENSTRAT
echo "par.BBF.EIGENSTRAT:"
cat par.BBF.EIGENSTRAT
echo
echo "------------------------------------------------------------------------"
echo

convertf -p par.BBF.EIGENSTRAT

if [[ $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo

# run EigenStrat

CMD="smartpca.perl -i $OutNam.eigenstratgeno -a $OutNam.snp -b $OutNam.ind -k 10 -o $OutNam.plus.HapMap.pca -p $OutNam.plus.HapMap.plot -e $OutNam.plus.HapMap.eval -l $OutNam.plus.HapMap.log -m 5 -t 2 -s 6.0"
echo $CMD
eval $CMD
if [ -f "$folder/$OutNam.plus.HapMap.eval" ] && [ -f "$folder/$OutNam.plus.HapMap.pca.evec"  ]; then
	#Rscript $HOME/CUMC/PCA_plot.R "${folder}/"  "$OutNam.plus.HapMap"
	python /home/local/users/jw/CUMC/Exome-Plots-Jiayao/PlotAncestryPCA.py -s $OutNam.plus.HapMap
fi

if [[ -e $BbfNam.fam.pcabkp ]]; then mv $BbfNam.fam.pcabkp $BbfNam.fam; fi
