import csv
import os
from optparse import OptionParser
import numpy as np
from collections import defaultdict
#from matplotlib.backends.backend_pdf import PdfPages   
import sys


#fname = 'PAH_dominamt_variants0716.csv'
#with open(DROPBOX + 'vqsr') as f:
#    vqsr = set(line.strip() for line in f.readlines())

#n = 324.0
parser = OptionParser()
parser.add_option("-c", "--csv", dest="csvfile",help="input csv file", metavar="csvfile")
parser.add_option("-n", "--numberofindividual", dest="Count",help="the numeber individuals in the csv file", metavar="Count")
#parser.add_option("-e", "--pop", dest="TXTfile",help="input sample list", metavar="TXTfile")
#parser.add_option("-x", "--exclude", dest="ExTXTfile",help="input excluded sample list", metavar="ExTXTfile")
#parser.add_option("-o", "--output", dest="OUTPUTfile",help="output file", metavar="OUTPUTfile")



(options, args) = parser.parse_args()
fname = options.csvfile
n=int(options.Count)
outname = fname.split('.csv')[0]+'2_filter.csv'
output_gene_count = fname.split('.csv')[0]+'_gene_counts2.csv'
output_variant_count = fname.split('.csv')[0]+'_variants_counts2.csv'
vqsr=set();


def check_rare_exon_damage(INFOdict, proband_gt, threshold = 0.001):
    '''find rare damage variant with population freq

    '''
    CodingCodes={'splicing', 'exonic', 'exonic,splicing'}
    LGD = {'splicing', 'frameshiftsubstitution', 'frameshift_insertion', 'frameshift_deletion', 'stopgain', 'stoploss'}    
    

    ExAcscore= max(map(float,[0]+[ rate  for rate in INFOdict.get('ExAC_ALL',str(threshold)).split(',') if rate not in {'NA','.'} ]))
    AF = max(map(float,[ rate  for rate in INFOdict.get('AF',str(threshold)).split(',') if rate != '.' ]))
    MutationFunct=str(INFOdict.get('Func.refGene','none').split(',')[0])
    MutationCalss = str(INFOdict.get('ExonicFunc.refGene','none').split(',')[0])
    MetaSVMprd = str(INFOdict.get('MetaSVM_pred','none').split(',')[0])
    if ',' in MutationCalss:
        gt = int(proband_gt.split('/')[-1])-1
        MutationCalss = MutationCalss.split(',')[gt]
                    
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







gene_counts = defaultdict(lambda: {'LGD':np.array([0, 0, 0]), 'DMIS':np.array([0, 0, 0])})
variant_counts = defaultdict(lambda: {'LGD':np.array([0, 0, 0]), 'DMIS':np.array([0, 0, 0])})

with open(fname, 'rU') as f:
    
    r = csv.reader(f)
    head = r.next()
    
    fw = open(outname, 'wb')
    w = csv.writer(fw)
    w.writerow(head)
    i = 0
    
    for line in r:
        variant = dict(zip(head, line))
        called =  int( variant['case_carrier']) +  int( variant['case_noncarrier'])
        af = sum(map( float, variant['AC'].split(',')))
        if variant['FILTER'] in vqsr: continue

        threshold = 0
        ExAcscore= max(map(float,[0]+[ rate  for rate in variant.get('ExAC_ALL',str(threshold)).split(',') if rate not in {'NA','.'} ]))
        if ExAcscore > 5*10**-5:
            continue

        
        
        if called > n * 0.15 and af < 10:
            w.writerow(line)
    
            homo, het, wild = 0, 0, 0
            gene, chrom, pos, ref, alt = variant['Gene.refGene'], variant['Gene.refGene'], variant['POS'] , variant['REF'], variant['ALT']
            proband_gt = variant['proband'].split('(')[1].split(':')[0]

            variant_type =  check_rare_exon_damage(variant, proband_gt, threshold = 5 * 10 ** -5 )
            variant_id = '_'.join([gene, chrom, pos, ref, alt])
            
            
            if proband_gt == '0/0':
                wild += 1
            elif proband_gt[0] == '0' and proband_gt[-1] != '0':
                het += 1
            elif proband_gt[0] == proband_gt[-1]:
                homo += 1

            if variant_type is 'LGD':
                    gene_counts[gene]['LGD'] = gene_counts[gene]['LGD'] + np.array([homo,het,wild])
                    variant_counts[variant_id]['LGD'] = variant_counts[variant_id]['LGD'] + np.array([homo,het,wild])

            elif variant_type is 'DMIS':
                gene_counts[gene]['DMIS'] = gene_counts[gene]['DMIS'] + np.array([homo,het,wild])
                variant_counts[variant_id]['DMIS'] = variant_counts[variant_id]['DMIS'] + np.array([homo,het,wild])

            
        
    fw.close()

fw = open(output_gene_count,'wb')
w = csv.writer(fw)
w.writerow(['Gene','LGD_homo', 'LGD_het','LGD_wt', 'DMIS_homo', 'DMIS_het', 'DMIS_wt','lgd_total', 'dmis_total', 'total'])
for gene, counts in gene_counts.items():
    lgd_counts = list(counts['LGD'])
    
    dmis_counts = list(counts['DMIS'])
    lgd_total = sum(lgd_counts[:2] )
    dmis_total= sum(dmis_counts[:2] )
    total = lgd_total + dmis_total
    w.writerow([gene] + lgd_counts + dmis_counts + [lgd_total, dmis_total, total] )
fw.close()

fw = open(output_variant_count ,'wb')
w = csv.writer(fw)
w.writerow(['Gene','LGD_homo', 'LGD_het','LGD_wt', 'DMIS_homo', 'DMIS_het', 'DMIS_wt'])
for gene, counts in variant_counts.items():
    w.writerow([gene] +list(counts['LGD']) + list(counts['DMIS']) )
fw.close()

