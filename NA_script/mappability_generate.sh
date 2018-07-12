#gem-indexer -T 10 -c dna -i $REF -o hg19_index
#gem-mappability -T 10 -I hg19_index.gem -l $1 -o "hg19_"$1"bp_mappability"  & #use screen for this step; take a very long time
#gem-2-wig -I hg19_index.gem -i "hg19_"$1"bp_mappability.mappability" -o "hg19_"$1"bp_mappability"


function wig2bed {
   fwig=$1
   fbed=$2
   echo "$fbed is generating!"
   awk 'BEGIN{
           OFS="\t";
           span=1;
           chr=0;
   }{
           if($1 ~/^var/){
                   split($2,a,"=");
                   chr=a[2];
                   if($4 ~/span/) {
                           split($4,a,"=");
                           span=a[2];
                   }else{
                           span=1
                   }
                   next
           }
            print chr,($1-1),($1+span-1),$2
         }' $fwig > $fbed
}

#### calculate the mappability based on the specified length of reads
###  it take very long time
function do_mappability {
        len=$1  ## set the length of the reads : 125
		ref=$2
		prefix=$3
        curpath=$(pwd)

        suff="${prefix}_${len}bp_mappability"

        cd $dst_map
        if [ ! -f "${suff}.wig" ]; then
				echo "generating $suff.wig!"
                gem-indexer -T 10 -c dna -i $ref -o ${prefix}_index
                gem-mappability -T 10 -I ${prefix}_index.gem -l $len -o $suff
                gem-2-wig -I ${prefix}_index.gem  -i ${suff}.mappability -o $suff
        fi
        if [ ! -f ${suff}.bed.gz ]; then
            wig2bed ${suff}.wig ${suff}.bed
                bgzip ${suff}.bed
                tabix -p bed ${suff}.bed.gz
        fi
        cd $curpath

}


dst_map="/home/local/ARCS/nz2274/Resources/mappability/"
echo "bash mappability_generate.sh 152 REF prefix_name"
do_mappability $1 $2 $3
