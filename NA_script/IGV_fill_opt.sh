pip_IGV="~/Pipeline/IGV/Generate_IGV_plots.sh"
igvscript="/home/local/ARCS/nz2274/Pipeline/NA_script/igvscript.sh"
len_read=76;
source /home/local/ARCS/nz2274/Pipeline/Exome_pipeline_scripts_GATKv3/WES_Pipeline_References.b37.sh

function bamout {

	fbamout=$1
    bam=$2
	ind=$3
	vars=$4
	interval="$ind.bed"
	
	awk -v id=$ind 'BEGIN{FS="\t";OFS="\t"}{if($1==id){print $2,$3-$len_read,$3+$len_read}}' $vars >$interval
	bamout="$fbamout/${ind}_bamout.bam" ## sandbox/NA12878_wgs_20.HC_out.bam 
	ovcf="$fbamout/${ind}_bamout.vcf" # "sandbox/NA12878_wgs_20_HC_calls_debug.vcf"
	
    if [ ! -f $bamout ]; then
	cmd="/home/local/users/jw/local/java/jre1.8.0_101/bin/java -jar $GATKJAR -T HaplotypeCaller \
	    -R $REF \
		-I $bam \
		-o $ovcf \
	    -bamout $bamout \
	    -forceActive -disableOptimizations \
	    -L  $interval"
   
         $cmd 
   fi
	#echo $bamout
} 

function Builtvariants {
	f_variants=$1; #"variants.txt" ### output
	f_ped=$2 #"../../Result/PAH_63.corrected.ped" ### ped file required
	f_var=$3 #"variant_info.txt"  ### ID\tCHR\tPOS\t
	f_bam=$4 #"../../BAM/bam.map.txt" ### ID\tBAM
	fbamout=$5
	while read -r line; do
		if [[ $line == "#"*  ]];then
			continue;     
		fi
		id=$(echo $line|awk '{print $1}')
		chr=$(echo $line|awk '{print $2}')
		lpos=$(echo $line|awk '{print $3}')
		FM=$(awk -v id=$id '{
			if($2==id){
				str="";
				if(length($3)!=1){
					print $3
				}
				if(length($4)!=1){
					print $4
			    }
			#print str
		}}' $f_ped)
		
		#########

		bams=""; #$(echo $bam_out|awk 'BEGIN{FS="/"}{print $(NF)}') 
		for s in $id $FM
		do
			n_bam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){n=split($2,b,"/");print b[n]}}')  
			bams=$n_bam","$bams
			#fullbam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){print $2}}'|head -n 1)
			#bam_out=$(bamout $fbamout $fullbam $chr $lpos $s)
	   	 	bams=$(echo $bam_out|awk 'BEGIN{FS="/"}{print $(NF)}')","$bams 
			
		done
		
		echo -e $chr"\t$lpos\t$bams" |sed 's/,$//g'  >>$f_variants
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
	f_ped=$3 #"../../Result/PAH_63.corrected.ped" ### ped file required
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
	
	for id in $(cut -f 1 $f_var|sort -k1,1d|uniq|grep -v '^#')
	do
		FM=$(awk -v id=$id '{
			if($2==id){
				str="";
				if(length($3)!=1){
					print $3
				}
				if(length($4)!=1){
					print $4
			    }
			#print str
		}}' $f_ped)
		bams=""; 
		for s in $id $FM
		do
			if [ $s != "" ];then 
				n_bam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){n=split($2,b,"/");print b[n]}}')  
				bams=$n_bam","$bams
				fullbam=$(egrep "^$s" $f_bam|awk -v id=$s '{if($1==id){print $2}}'|head -n 1)
				bamout $fbamout $fullbam $s $f_var
		   	 	bam_out="$fbamout/${ind}_bamout.bam" 
				bams=$(echo $bam_out|awk 'BEGIN{FS="/"}{print $(NF)}')","$bams 
			    echo $bam_out >> $fbamlist;
			fi
		done
		echo "BAM\t"$bams
	   awk -v id=$id -v bams=$bams 'BEGIN{FS="\t";OFS="\t"}{if($1==id){print $2,$3,bams}}' $f_var  >> $fvariants
	   
	done
	

	sed -i 's/,$//g' $fvariants

#	buildbam $fbamout $fbamlist
	cut -f 2 $f_bam >> $fbamlist
    #buildbam $fbam $fbamlist
	${igvscript}  -v $fvariants -b $fbamlist
	echo "${pip_IGV} -v $fvariants -b $fbamlist"

}

f_ped=$1; # example: "PAH_Vanderbilt/Result/PAH_63.corrected.ped" ### ped file required
f_var=$2; # example: "variant2_info.txt"  ### ID\tCHR\tPOS\t
f_bam=$3; # example: "PAH_Vanderbilt/BAM/bam.map.txt" ### ID\tBAM

if [ ! -e $f_ped  ]; then
	echo $fped" cannot find!";
	exit;
fi

if [  ${f_ped: -4} == ".ped"  ]
then
   #echo "$f_ped. the order of input is: ped variant.list bam.list "
   #exit
#fi
#exit;
call_igv $f_var $f_bam $f_ped 

else 
  echo "input ped variant bam"
fi
### bash IGV_fill.sh ped var bam
