vcf=$1; #  "input.vcf"
bed=$2; #"bed_file_describing_the_range.bed"
output=$3; #"output_prefix"
vcftools --vcf  $vcf --bed $bed --out $output --recode --keep-INFO-all