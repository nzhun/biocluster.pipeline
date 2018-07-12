vcf=$1
fname=$(basename $vcf)
fbname=${fname%.*}
nohup zcat  $vcf |python ~/Pipeline/scripts/lift_over.py  --format vcf > $fbname.hg38.vcf  2> nohup.$fbname.log &
