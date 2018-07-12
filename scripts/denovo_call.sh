echo "input: vcf ped name"
vcf=$1
ped=$2
outNm=$3
folder=$(dirname $ped)
python ~/Pipeline/scripts/adhoc.1.vcf2trio.py -v $vcf  -p $ped
python ~/Pipeline/scripts/adhoc.2.trio2denovo.py -v $folder/vcf_trio/
python ~/Pipeline/scripts/adhoc.3.variants2csv.py -v $folder/vcf_rare_denovo/ -o $outNm
python ~/Pipeline/scripts/adhoc.4.denovo_filter.py -c $outNm  -f  ~/Pipeline/scripts/denovo_soft2.yml
