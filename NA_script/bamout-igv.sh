function bamout {

	fbamout=$1
    bam=$2
	chr=$3
	pos=$4
	ind=$5
	beg=$((pos - len_read))
	end=$((pos + len_read))
	bamout="$fbamout/${ind}_${chr}_${pos}_bamout.bam" ## sandbox/NA12878_wgs_20.HC_out.bam 
	ovcf="$fbamout/${ind}_${chr}_${pos}_bamout.vcf" # "sandbox/NA12878_wgs_20_HC_calls_debug.vcf"
	interval="${chr}:${beg}-${end}"
#	echo "a"$GATKJAR
    if [ ! -f $bamout ]; then
	java -jar $GATKJAR -T HaplotypeCaller \
	    -R $REF \
		-I $bam \
		-o $ovcf \
	    -bamout $bamout \
	    -forceActive -disableOptimizations \
	    -L  $interval 
   fi
	echo $bamout
} 

fvar="variants.txt"
