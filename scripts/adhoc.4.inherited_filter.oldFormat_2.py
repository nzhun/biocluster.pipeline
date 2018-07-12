
from __future__ import division
import csv
import os
from optparse import OptionParser
import yaml


def reads_filter(format, proband, parents, parameter): # R,A used for checking indel AC is [1,1]
    reads_pass = True
    ID, geno =  proband.split('(')
    father,mother = parents.split('),')
    fatherID,fathergeno = father.split('(')
    motherID,mothergeno = mother[:-1].split('(')

    
    proband_format = dict(zip(format.split(':'), geno[:-1].split(':')))
    GT,AD,DP,GQ,PL =  proband_format['GT'], proband_format['AD'], proband_format['DP'], proband_format['GQ'], proband_format['PL']
    
        
    ref_idx, alt_idx = map(int,GT.split('/'))
    
    if  DP in {'.','0'}:
        DP = 1  

    # prband filter
    proband_AD = map(int,AD.split(','))
    proband_DP = float(DP)
    proband_PL = int(PL.split(',')[0])
    proband_maf = proband_AD[alt_idx]/proband_DP

    if proband_PL < parameter['reads_filter']['min_proband_PL']:
        reads_pass = False
    if proband_AD[alt_idx] <  parameter['reads_filter']['min_proband_AD']:
        reads_pass = False

    if proband_AD[alt_idx] >= 10:
        if proband_maf < parameter['reads_filter']['min_proband_alt_freq_tier2']:
            reads_pass = False
    else:
        if proband_maf < parameter['reads_filter']['min_proband_alt_freq_tier1']:
            reads_pass = False
        
    return  reads_pass
        
def info_filter(variants, snps, parameter):
    '''
    functions that flag the variant of FS, QD
    '''
    p = True

    # freq filter
    ExACfreq = max(map(float,[0] + [rate for rate in variants.get('ExACfreq','0').split(',') if rate not in {'NA','.'} ]))
    ESPfreq = max(map(float,[0] + [rate for rate in variants.get('ESPfreq','0').split(',') if rate not in {'NA','.'} ]))
    KGsfreq = max(map(float,[0] + [rate for rate in variants.get('1KGfreq','0').split(',') if rate not in {'NA','.'} ]))
    if ExACfreq > parameter['gene_filter']['max_population_AF']['ExAC'] or \
       ESPfreq > parameter['gene_filter']['max_population_AF']['ESP'] or \
       KGsfreq > parameter['gene_filter']['max_population_AF']['1KG']:
       p = False

    # exon filter
    if parameter['gene_filter']['exon']:
        if variants['VarFunc'] not in parameter['gene_filter']['exon_flag']:
            p = False

    if variants['QD'] == 'NA':
        variants['QD'] = 100
    if snps:
        if float(variants['FS']) > parameter['snps']['max_FS']:
            p = False
        if variants['QD'] == 'NA' or float(variants['QD']) < parameter['snps']['min_QD']:
            p = False 
    else:
        if float(variants['FS']) > parameter['indels']['max_FS']:
            p = False
        if float(variants['QD']) < parameter['indels']['min_QD']:
            p = False   
        if variants['ReadPosRankSum'] != 'NA' and float(variants['ReadPosRankSum']) < parameter['indels']['min_ReadPosRankSum']:
            p = False         
    return p
 

def variants_filter(filename, parameter):
    
    with open(filename,'Ur') as denovo_f:
        denovo_r = csv.reader(denovo_f)
        head = denovo_r.next()
        folder=os.path.dirname(filename)
        prefix=os.path.basename(filename)
        denovo_snp, denovo_indel =[], []
       # prefix='.'.join(filename.rsplit('.', 1)[:-1] )
        print "filtered result writes to "+ folder+"/"+prefix  + '_filtered.csv'
        denovo_write= open(folder+'/'+prefix + '_filtered.csv','wb')
        w= csv.writer(denovo_write)
        w.writerow(head) 
        
        idx = dict(zip(head, range(len(head))))
        for line in denovo_r:


            variant_pass = True
            
            
            variant = dict(zip(head, line))
            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']

            
            format = variant['FORMAT']
            probands = variant['proband'].split(';')
            parents = variant['parents'].split(';')
            AC = variant['AC'].split(',') 
            AF = max(map(float, variant['AF'].split(',')))
           # print str(AF)+"\t"+str(AC)+"\n"
            if 'max_cohort_AF' in parameter['gene_filter']:
                if AF > parameter['gene_filter']['max_cohort_AF'] : 
                    variant_pass = False

            
            # remove multi-allel sites
            if len(alt.split(',')) > parameter['gene_filter']['max_Multiallelic'] : 
                variant_pass = False
            #print str(len(alt.split(',')))
            # remove know problem genes
            if chrom in  parameter['gene_filter']['excluded_chrom']:
                variant_pass = False
            for e in parameter['gene_filter']['excluded_gene']:
                if e in variant['GeneName']:
                    variant_pass = False
                    break
            # remove high segdup score  
        
            segdup = 0 if variant.get('genomicSuperDups') in {'NA','.'} else float(variant.get('genomicSuperDups',':0').split(':')[1].split('-')[0])
            #print segdup
            if segdup > parameter['gene_filter']['max_seqdup']: 
                variant_pass = False
            
            for i in range(len(probands)) :
                newline = list(line)
                proband = probands[i]
                newline[idx['proband']] = proband
                parent = parents[i] 
                newline[idx['parents']] = parent
                ID, geno = proband.split('(')
                GT = geno[:-1].split(':')[0]
                _, alt_idx = map(int,GT.split('/'))
                alt_idx=1;
               # print alt
                newline[idx['ALT']] = alt.split(',')[alt_idx-1]                
                real_AC = int(AC[alt_idx-1])
                if 'max_cohort_AC' in parameter['gene_filter'] :
                    if real_AC > parameter['gene_filter']['max_cohort_AC']:
                        variant_pass = False

                snp = len(newline[idx['ALT']]) == len(ref) and newline[idx['ALT']] != '*'

               # print "1  "+str(real_AC)
                
                if not info_filter(variant, snp, parameter):
                    variant_pass = False
              #  print "2  "+str(variant_pass)
                variant_id = '_'.join([ID,chrom,pos])
                if len(variant['VarClass'].split(',')) > 1:
                    newline[idx['VarClass']] = variant['VarClass'].split(',')[alt_idx-1]
              #  print "3  "+str(variant_pass)
                if snp: 
                    if variant_pass and reads_filter(format, proband, parent, parameter) :
                        denovo_snp.append(variant_id)
                        w.writerow(newline)   
                else:
                    if variant_pass and reads_filter(format, proband, parent, parameter):
                        denovo_indel.append( variant_id)
                        w.writerow(newline)               

        denovo_write.close()
        denovo_f.close()
    return denovo_snp, denovo_indel

 
parser = OptionParser()
parser.add_option("-c", "--csv", dest="CSVfile",help="input CSV file", metavar="CSVfile")
parser.add_option("-f", "--cfg", dest="CFGfile",help="input de novo filter configure file", metavar="CFGfile")
(options, args) = parser.parse_args()
csv_name = options.CSVfile
cfg_name = options.CFGfile


csv_name  = os.path.abspath(csv_name)

with open(cfg_name, 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

print '-' * 50
for section in cfg:
    print(section)
    for key, value in cfg[section].items():
        print '     %s: %s' % (key, value)
print '-' * 50

denovo_snp, denovo_indel = variants_filter(csv_name,cfg)
print '%s denovo snp, %s denovo indels' % (len(denovo_snp), len(denovo_indel))


