
from __future__ import division
import os
from optparse import OptionParser
import gzip
import csv
from collections import defaultdict
import numpy as np
import tabix

parser = OptionParser()
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-p", "--ped", dest="PEDfile",help="input PED file", metavar="PEDfile")
parser.add_option("-e", "--pop", dest="TXTfile",help="input sample list", metavar="TXTfile")
parser.add_option("-x", "--exclude", dest="ExTXTfile",help="input excluded sample list", metavar="ExTXTfile")
parser.add_option("-o", "--output", dest="OUTPUTfile",help="output file", metavar="OUTPUTfile")



(options, args) = parser.parse_args()
vcf_name = options.VCFfile
ped_name = options.PEDfile
european = options.TXTfile
exclude =  options.ExTXTfile
output =  options.OUTPUTfile

#PROJECT = '/home/local/ARCS/hq2130/WES/PAH/'
#vcf_name = PROJECT + 'vcf0716/PAH_control_VQSR.vcf'
#european = PROJECT + 'src/european.txt'
#ped_name = PROJECT + 'src/PAH.ped'
#exclude = PROJECT + 'src/excluded_sample.txt'

output_variants = output
output_gene_count =  output+'_gene_counts.csv'
output_variant_count = output+'_variants_counts.csv'

def check_rare_exon_damage(INFOstring, threshold = 0.001):
    '''find rare damage variant with population freq

    '''
    INFOcolumnList=INFOstring.split(";")
    INFOdict={}
    for element in INFOcolumnList:
        if '=' in element:
            FieldName,FieldValue=element.split('=',1)
            INFOdict[FieldName]=FieldValue

    threshold = threshold
    CodingCodes={'splicing', 'exonic', 'exonic,splicing'}
    LGD = {'splicing', 'frameshiftsubstitution', 'frameshift_insertion', 'frameshift_deletion', 'stopgain', 'stoploss'}    
    
   # ExACfreq = max(map(float,[0] + [rate for rate in variants.get('ExAC_ALL','0').split(',') if rate not in {'NA','.'} ]))
#    ESPfreq = max(map(float,[0] + [rate for rate in variants.get('esp6500siv2_all','0').split(',') if rate not in {'NA','.'} ]))
 #   KGsfreq = max(map(float,[0] + [rate for rate in variants.get('1000g2015aug_all','0').split(',') if rate not in {'NA','.'} ]))
    
    ExAcscore= max(map(float,[0]+[ rate  for rate in INFOdict.get('ExAC_ALL',str(threshold)).split(',') if rate not in {'NA','.'} ]))
    AF = max(map(float,[ rate  for rate in INFOdict.get('AF',str(threshold)).split(',') if rate != '.' ]))
    MutationFunct=str(INFOdict.get('Func.refGene','none').split(',')[0])
    MutationCalss = str(INFOdict.get('ExonicFunc.refGene','none').split(',')[0])
    MetaSVMprd = str(INFOdict.get('MetaSVM_pred','none').split(',')[0])
    LGD = (MutationCalss in LGD) or ( MutationFunct == 'splicing')
    Dmis = (MutationCalss == 'nonsynonymous_SNV' and MetaSVMprd == 'D')
    gene = INFOdict.get('GeneName','NA')

    
    if ExAcscore <= threshold and MutationFunct in CodingCodes:

        if LGD is True:
            return 'LGD'
        elif Dmis is True:
            return 'DMIS'
        else:
            return 'not damage'
    else:
        return 'fail'







vqsr_exclude = set() ## [ "VQSRTrancheINDEL99.90to100.00","VQSRTrancheINDEL99.90to100.00+","VQSRTrancheSNP99.90to100.00","VQSRTrancheSNP99.90to100.00+"]

# pick european samples
with open(ped_name) as f:
    cases = set(line.strip().split()[1] for line in f.readlines())
with open(european) as f:
    euro = set(line.strip() for line in f.readlines())
with open(exclude) as f:
    exclude = set(line.strip() for line in f.readlines())
subset = cases & euro - exclude

with gzip.open(vcf_name, 'r') as f:
    for line in f: 
        if line[:6]=='#CHROM': 
            samples = line.split()
            break

# pick index of proband 
proband_index=[]
for proband in subset:
    if proband in samples:
        proband_index.append(samples.index(proband))
print len(proband_index) , 'european total in domianmt model analysis'


# pedigree index for finding variants
sample_index={}
with open(ped_name) as f:
    for line in f:
        family, proband, father, mother = line.split()[:4]
        if proband in subset:
            if proband in samples and father in samples and mother in samples:
                sample_index[proband]=[father,mother,samples.index(proband),samples.index(father),samples.index(mother)]
            elif proband in samples and father in samples:
                sample_index[proband]=[father,'unknown',samples.index(proband),samples.index(father),samples.index('ID')]
            elif proband in samples and mother in samples:
                sample_index[proband]=['unknown', mother, samples.index(proband),samples.index('ID'), samples.index(mother)]
            elif proband in samples:
                sample_index[proband]=['unknown','unknown',samples.index(proband),samples.index('ID'),samples.index('ID')]     
                

    
gene_counts = defaultdict(lambda: {'LGD':np.array([0, 0, 0]), 'DMIS':np.array([0, 0, 0])})
variant_counts = defaultdict(lambda: {'LGD':np.array([0, 0, 0]), 'DMIS':np.array([0, 0, 0])})



print 'Start filtering variants '
if vcf_name.endswith(".gz") :
    f=gzip.open(vcf_name)
else:
    f=open(vcf_name)
#with gzip.open(vcf_name)  as f:
    
fw = open(output_variants ,'wb')
w = csv.writer(fw)
head_info = ['Gene.refGene','Func.refGene', 'ExonicFunc.refGene','ABHet', 'ABHom', 'AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'DB', 'DP', 'DS', 'Dels', 'END', 'FS', 'GQ_MEAN', \
'GQ_STDDEV', 'HRun', 'HapMapV3', 'HaplotypeScore', 'InbreedingCoeff', 'MLEAC', 'MLEAF', 'MQ', 'MQ0', 'MQRankSum', 'NCC', \
'OND', 'QD', 'ReadPosRankSum', 'SOR', 'AAChange.refGene', 'esp6500siv2_all', 'esp6500siv2_aa', 'esp6500siv2_ea', '1000g2015aug_all', '1000g2015aug_eur', \
'1000g2015aug_amr', '1000g2015aug_eas', '1000g2015aug_afr', '1000g2015aug_sas', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', \
'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score',\
'Polyphen2_HVAR_pred', 'MutationTaster_score', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_pred', 'MetaSVM_score', 'MetaSVM_pred', 'CADD13_RawScore', 'CADD13_PHRED', 'GERP++_RS', 'phyloP20way_mammalian', 'SiPhy_29way_logOdds', \
'FATHMM_coding', 'FATHMM_noncoding', 'GWAVA_region_score', 'GWAVA_tss_score', 'GWAVA_unmatched_score', \
'cosmic70' , 'genomicSuperDups','MCAP']
w.writerow(['case_carrier', 'case_noncarrier','proband', 'parents', 'FORMAT', 'CHROM', 'POS', 'ID',  'REF', 'ALT', 'QUAL','FILTER'] + head_info)


for line in f:
    if line[0] == '#':
        continue
    else:
        data =  line.split('\t')
        vcf_size=len(data)-9
        INFOstring = data[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                if FieldName !="VQSR" and FieldName in INFOdict :
                    FieldValue=INFOdict.get(FieldName)+","+FieldValue
                INFOdict[FieldName]=FieldValue
        FORMAT = data[8]
        if int(INFOdict.get('AN')) < vcf_size*0.9*2: continue
        if int(INFOdict.get('DP')) < 15*vcf_size : continue 
        gene = INFOdict.get('Gene.refGene','NA')
        AC = INFOdict.get('AC').split(',')
        allows=[];
        for i in range(0,len(AC)):
            if int(AC[i])<vcf_size*0.05 :
                allows.append(str(i+1));
        if len(allows)==0 : continue
        if len(AC) > 2: continue
       # if min(map(int,AC)) > int(INFOdict.get('AN'))*0.05: continue
       # if max(map(int,AC)) > 15: continue
        if 'MUC' in gene: continue
        if 'HLA' in gene: continue
        if data[6] in vqsr_exclude: continue
        
        CodingCodes={'splicing', 'exonic', 'exonic,splicing'}
        MutationFunct=str(INFOdict.get('Func.refGene','none').split(',')[0])
        if MutationFunct not in CodingCodes:
            continue
      #  print  INFOdict.get('genomicSuperDups')     
        segdup = 0 if INFOdict.get('genomicSuperDups') in {'NA','.','.,.'} else float(INFOdict.get('genomicSuperDups',':0').split(':')[1].split('-')[0])
        if segdup > 0.95: 
            continue

        # counts carrier number for this variant
        homo, het, wild = 0, 0, 0
        variant_type =  check_rare_exon_damage(INFOstring, threshold = 5 * 10 ** -5)
        variant_id = '_'.join([gene, data[0], data[1], data[3], data[4]])


        if True: # keep all type of variants variant_type in {'LGD','DMIS'}:
            
            for index in proband_index:
                info = dict(zip(data[8].split(':'),data[index].split(':')))
                if 'GQ' not in info: continue
                if 'DP' not in info: continue
                if info['DP'] == '.': continue
                if info['GQ'] == '.': continue
                ads=info['AD'].split(",");
                maxad=int(ads[0])
                for l in range(1,len(ads)) :
                    if int(ads[l]) > maxad :
                        maxad=int(ads[l])
                if int(info['GQ']) >= 60 and int(info['DP']) >= 8:
                    if info['GT'] == '0/0':
                        wild += 1

                    elif info['GT'][0] == '0' and info['GT'][-1] != '0' :
                        if info['GT'][-1] in allows and maxad < (int(info['DP']) -3) :
                            het += 1

                    elif info['GT'][0] == info['GT'][-1]:
                        if info['GT'][-1] in allows :
                            homo += 1
                        
        if homo + het > 0:
            
            if variant_type is 'LGD':
                gene_counts[gene]['LGD'] = gene_counts[gene]['LGD'] + np.array([homo,het,wild])
                variant_counts[variant_id]['LGD'] = variant_counts[variant_id]['LGD'] + np.array([homo,het,wild])

            if variant_type is 'DMIS':
                gene_counts[gene]['DMIS'] = gene_counts[gene]['DMIS'] + np.array([homo,het,wild])
                variant_counts[variant_id]['DMIS'] = variant_counts[variant_id]['DMIS'] + np.array([homo,het,wild])
                
            # write varaint details into csv
            new_info = []
            for e in head_info:
                new_info.append(str(INFOdict.get(e,'NA')))

            for proband in sample_index:
                father, mother, proband_index2, father_index, mother_index = sample_index[proband]
                proband_gt =  data[proband_index2].split(':')[0]
                iac = max(map(int,[ allele  for allele in proband_gt.split('/') if allele not in {'.'} ])+[0])
                if str(iac) in allows and proband_gt not in {'0/0', './.'}:
                    info = dict(zip(FORMAT.split(':'),data[proband_index2].split(':')))
                    if 'GQ' not in info: continue
                    if 'DP' not in info: continue
                    if int(info['GQ']) >= 60 and int(info['DP']) >= 8:
                    
                        proband_info = proband + '(' + data[proband_index2] + ')'
                        parents_info = father + '(' + data[father_index] + ');' + mother + '(' + data[mother_index] + ')'
                        w.writerow( [het + homo, wild, proband_info, parents_info, FORMAT] + data[:7] + new_info)
fw.close()
                
                


fw = open(output_gene_count,'wb')
w = csv.writer(fw)
w.writerow(['Gene','LGD_homo', 'LGD_het','LGD_wt', 'DMIS_homo', 'DMIS_het', 'DMIS_wt'])
for gene, counts in gene_counts.items():
    w.writerow([gene] +list(counts['LGD']) + list(counts['DMIS']) )
fw.close()

fw = open(output_variant_count ,'wb')
w = csv.writer(fw)
w.writerow(['Gene','LGD_homo', 'LGD_het','LGD_wt', 'DMIS_homo', 'DMIS_het', 'DMIS_wt'])
for gene, counts in variant_counts.items():
    w.writerow([gene] +list(counts['LGD']) + list(counts['DMIS']) )
fw.close()
