# /home/nz2274/UK/pah-icon/vcf/gvcf.list  /home/nz2274/UK/pah-icon/src/vcr_xgen_union.hg19.bed 2000 
gvcf=$1
bed=$2
len=$3
total=$(wc $bed|awk '{print $1}')
echo $total
#seq 0 $((total+len+1))  $len
for ((i=0;i< $((total+len));i=i+len)); do 
echo $i
#qsub -N genotype.$i /home/nz2274/Pipeline/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh  -i gvcf.list  -r /home/nz2274/Pipeline/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh  -t /home/nz2274/UK/pah-icon/src/vcr_xgen_union.hg19.bed -a $i -s $len
qsub -q shenlab.q@compute-0-62.local -t 1-2 -N genotype.$i /home/nz2274/Pipeline/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh  -i  $gvcf  -r /home/nz2274/Pipeline/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh  -t $bed -a $i -s $len
done

