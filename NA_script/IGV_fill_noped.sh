pip_IGV="~/Pipeline/IGV/Generate_IGV_plots.sh"
len_read=125;
source /home/local/ARCS/nz2274/Pipeline/Exome_pipeline_scripts_GATKv3/WES_Pipeline_References.b37.sh

function bamout {

	fbamout=$1
    bam=$2
	chr=$3
	pos=$4
	ind=$5
	beg=$((pos - len_read - len_read ))
	end=$((pos + len_read + len_read ))
	bamout="$fbamout/${ind}_${chr}_${beg}_${end}_bamout.bam" ## sandbox/NA12878_wgs_20.HC_out.bam 
	ovcf="$fbamout/${ind}_${chr}_${beg}_${end}_bamout.vcf" # "sandbox/NA12878_wgs_20_HC_calls_debug.vcf"
	interval="${chr}:${beg}-${end}"
#	echo "a"$GATKJAR
    if [ ! -f $bamout ]; then
	/home/local/users/jw/local/java/jre1.8.0_101/bin/java -jar $GATKJAR -T HaplotypeCaller \
	    -R $REF \
		-I $bam \
		-o $ovcf \
	    -bamout $bamout \
	    -forceActive -disableOptimizations \
	    -L  $interval 
   fi
	echo $bamout
} 

function Builtvariants {
	f_variants=$1; #"variants.txt" ### output
	#f_ped=$2 #"../../Result/PAH_63.corrected.ped" ### ped file required
	f_var=$2 #"variant_info.txt"  ### ID\tCHR\tPOS\t
	f_bam=$3 #"../../BAM/bam.map.txt" ### ID\tBAM
	fbamout=$4
	while read -r line; do
		if [[ $line == "#"*  ]];then
			continue;     
		fi
		id=$(echo $line|awk '{print $1}')
		chr=$(echo $line|awk '{print $2}')
		pos=$(echo $line|awk '{print $3}')
		#FM=$(awk -v id=$id '{
		#	if($2==id){
	#			str="";
#				if($3!="-"&&$3!="x"){#
			#		print $3
		#		}
	#			if($4!="-" && $4!="x"){
#					print $4#
#			    }
			#print str
#		}}' $f_ped)
		
		#########

		bams=""; #$(echo $bam_out|awk 'BEGIN{FS="/"}{print $(NF)}') 
		for s in $id # $FM
		do
			n_bam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){n=split($2,b,"/");print b[n]}}')  
			bams=$n_bam","$bams
			fullbam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){print $2}}')
			bam_out=$(bamout $fbamout $fullbam $chr $pos $s)
	   	 	bams=$(echo $bam_out|awk 'BEGIN{FS="/"}{print $(NF)}')","$bams 
			
		done
		
		echo -e $chr"\t$pos\t$bams" |sed 's/,$//g'  >>$f_variants
	done < $f_var
    #echo $f_variants       
}


function buildbam {
	folder=$1;
	fout=$2
	for f in $folder/*.bam
	do
		readlink -f $f >>$fout
	done
}

function call_igv {
	#f_ped=$3 #"../../Result/PAH_63.corrected.ped" ### ped file required
	f_var=$1 #"variant_info.txt"  ### ID\tCHR\tPOS\t
	f_bam=$2 #"../../BAM/bam.map.txt" ### ID\tBAM
	echo "$GATKJAR"
	fbamlist=$(basename $f_var)"_bam.txt"
	if [ -e $fbamlist ]; then
	   rm $fbamlist
    fi
	fvariants=$(basename $f_var)"_variants.txt"
	if [ -e $fvariants ]; then 
		rm $fvariants
	fi
	fbamout="BAMOUT"
#	fbam="."
	mkdir $fbamout
	
	

	Builtvariants $fvariants  $f_var $f_bam $fbamout

	buildbam $fbamout $fbamlist
	cut -f 2 $f_bam >> $fbamlist
    #buildbam $fbam $fbamlist
	echo "${pip_IGV} -v $fvariants -b $fbamlist"
}

#f_ped=$1; # example: "PAH_Vanderbilt/Result/PAH_63.corrected.ped" ### ped file required
f_var=$1; # example: "variant2_info.txt"  ### ID\tCHR\tPOS\t
f_bam=$2; # example: "PAH_Vanderbilt/BAM/bam.map.txt" ### ID\tBAM

call_igv $f_var $f_bam  

