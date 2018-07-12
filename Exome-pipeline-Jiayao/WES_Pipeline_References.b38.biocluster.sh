## Resource Directories
FOLDS="/home/yufengshen"
export EXOMPPLN="$FOLDS/CUMC/Exome-pipeline-Jiayao/" # Directory containing pipeline shell scripts
export EXOMRES="$FOLDS/resources" # Directory containing resources/references for pipeline

#jar files  and directories for software
GATKJAR="$FOLDS/bin/GenomeAnalysisTK.jar" #Current GATK jar file
PICARD="$FOLDS/bin/picard.jar" #directory containing Picard jar files
SNPEFF="$FOLDS/bin/snpEff.jar" # Current snpEff jar file
Trimmomatic="$FOLDS//software_pkg/Trimmomatic-0.36/trimmomatic-0.36.jar" # Current snpEff jar file

## References
export BUILD="hg38" # shorthand for build
#export REF="$EXOMRES/reference_genomes/GRCh38/Homo_sapiens_assembly38.fasta" 
export REF="/share/data/resources/hg38/references/GRCh38.no_alt_analysis_set.fa" 
export HAPMAP="/share/data/resources/hg38/references/hapmap_3.3.hg38.vcf" # hapmap vcf from GATK
export INDEL="/share/data/resources/hg38/references/Mills_and_1000G_gold_standard.indels.hg38.vcf" # Gold standard INDEL reference from GATK
export TGVCF="/share/data/resources/hg38/references/1000G_omni2.5.hg38.vcf" 
export DBSNP="/share/data/resources/hg38/references/dbsnp_146.hg38.vcf" # dbSNP vcf from GATK
#export DBSNP="$EXOMRES/references/GRCh38/dbsnp_138.hg38.excluding_sites_after_129.sorted.cleaned.vcf " # dbSNP vcf from GATK
export ONEKG="/share/data/resources/hg38/references/1000G_phase1.snps.high_confidence.hg38.vcf" # 1000 genome SNPs vcf
export ANNOVAR="$FOLDS/software_pkg/annovar"
export ANNHDB="/share/data/resources/ANNOVAR_DATA/humandb" #Location of annovar databases

#Capture Kit Target Files
export TGTCODES="AgtV2:AgtV4:AgtV5:AgtV5UTR:NbgV2:NbgV3:IllTS:BigTgt:RefSeq:VCRv2"
export Adapter_TruSeq3_PE="/home/yufengshen/software_pkg/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
export Adapter_TruSeq3_SE="/home/yufengshen/software_pkg/Trimmomatic-0.36/adapters/TruSeq3-SE.fa"
