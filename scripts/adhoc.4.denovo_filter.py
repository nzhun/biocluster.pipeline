
from __future__ import division
import csv
import os
import re
from optparse import OptionParser
import yaml


def check_damage(variants, varclass):
    '''
    function return if a variant if damage 
    '''
    damage = False

    damage_key = {'frameshift','stop','splicing'}

    for key in damage_key:
        if varclass.startswith(key): 
            damage = True
    if 'nonsynonymousSNV' in  varclass:
        if 'D' in variants['MetaSVMprd']:
            damage = True
    return damage




def reads_filter(format, proband, parents, parameter): # R,A used for checking indel AC is [1,1]
    reads_pass = True
    ID, geno =  proband.split('(')
    father,mother = parents.split('),')
    fatherID,fathergeno = father.split('(')
    motherID,mothergeno = mother[:-1].split('(')
    #print parents
    
    proband_format = dict(zip(format.split(':'), geno[:-1].split(':')))
    GT,AD,DP,GQ,PL =  proband_format['GT'], proband_format['AD'], proband_format['DP'], proband_format['GQ'], proband_format['PL']
    father_format = dict(zip(format.split(':'), fathergeno.split(':')))
    fGT,fAD,fDP,fGQ,fPL =  father_format['GT'], father_format['AD'], father_format['DP'], father_format['GQ'], father_format['PL']
    mother_format = dict(zip(format.split(':'), mothergeno[:-1].split(':')))
    mGT,mAD,mDP,mGQ,mPL  =  mother_format['GT'], mother_format['AD'], mother_format['DP'], mother_format['GQ'], mother_format['PL']

        
    ref_idx, alt_idx = map(int,GT.split('/'))
    
    # parents filter
    if  mDP == '.':
       # print mother
        mDP = 1
    if  fDP == '.':
        fDP = 1
    if  DP in {'.','0'}:
        DP = 1  
                
    if float(mDP) < 1:
        mDP = 1
    if float(fDP) < 1:
        fDP = 1
    Fmaf = int(fAD.split(',')[alt_idx])/float(fDP)
    Mmaf = int(mAD.split(',')[alt_idx])/float(mDP) 
    # <<<<<<<<<<<<<<<<<<<<<<

     # parents filter
    if max(Fmaf, Mmaf) > parameter['reads_filter']['max_parents_alt_freq']:
        reads_pass = False
    if min(int(mDP), int(fDP)) < parameter['reads_filter']['min_parents_DP']:
        reads_pass = False
    if min(int(mGQ), int(fGQ)) < parameter['reads_filter']['min_parents_GQ']:
        reads_pass = False
   # print reads_pass   
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
   # print "QD: "+str(float(variants['QD']))+"\t"+"FS\t"+str(float(variants['FS']))+" ReadPos\t"+str(variants['ReadPosRankSum'])
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
    
    #print "p: "+ str(p)   
    return p
 

def variants_filter(filename, parameter):
    
    with open(filename,'Ur') as denovo_f:
        denovo_r = csv.reader(denovo_f)
        head = denovo_r.next()

        denovo_snp, denovo_indel =[], []
        prefix='.'.join(filename.rsplit('.', 1)[:-1] )
        denovo_write= open(prefix+ '_filtered.csv','wb')
        w= csv.writer(denovo_write)
        w.writerow(head) 

        denovo2_write= open(prefix + '_lessconfident_damage_filtered.csv','wb')
        w2 = csv.writer(denovo2_write)
        w2.writerow(head) 

        
        idx = dict(zip(head, range(len(head))))
        for line in denovo_r:


            variant_pass = True
            
            
            variant = dict(zip(head, line))
            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']
            alts=alt.split(',')
            
            format = variant['FORMAT']
            probands = variant['proband'].split(';')
            parents = variant['parents'].split(';')
            AC = variant['AC'].split(',') 
            
            # remove multi-allel sites
            if len(alt.split(',')) > parameter['gene_filter']['max_Multiallelic'] : 
                variant_pass = False
            
            # remove know problem genes
            if chrom in  parameter['gene_filter']['excluded_chrom']:
                variant_pass = False
            for e in parameter['gene_filter']['excluded_gene']:
                if e in variant['Gene.refGene']:
                    variant_pass = False
                    break

            # remove high segdup score  
            #print variant.get('genomicSuperDups')
            segdup = 0 if variant.get('genomicSuperDups') in {'NA','.'} or re.match('.,',variant.get('genomicSuperDups'))  else float(variant.get('genomicSuperDups',':0').split(':')[1].split(',')[0].split('-')[0]) #
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
                alt_idx=1
                newline[idx['ALT']] = alts[alt_idx-1]   
                real_AC = int(AC[alt_idx-1])
                if real_AC > parameter['gene_filter']['max_cohort_AC']:
                    variant_pass = False
       
                total_AC = sum(map(int, AC))
                if total_AC > 2 * parameter['gene_filter']['max_cohort_AC']:
                    variant_pass = False
            
                snp = len(newline[idx['ALT']]) == len(ref) and newline[idx['ALT']] != '*'

               # if real_AC > parameter ['snps']['max_AC']:
                #    variant_pass=False
                #    print 'AC: %d' %(real_AC)
                if not info_filter(variant, snp, parameter):
                    variant_pass = False
             
                variant_id = '_'.join([ID,chrom,pos])
                
                # change varclass format
                varclass = variant['ExonicFunc.refGene']
               
                if len(variant['ExonicFunc.refGene'].split(',')) > 1:
                    varclass = variant['ExonicFunc.refGene'].split(',')[alt_idx-1]
                    newline[idx['ExonicFunc.refGene']] = varclass
                if snp: 
                    if variant_pass and reads_filter(format, proband, parent, parameter):
                        denovo_snp.append(variant_id)
                        w.writerow(newline)  
                    elif variant_pass and check_damage(variant, varclass):
                        w2.writerow(newline)  

                else:
                    if variant_pass and reads_filter(format, proband, parent, parameter):
                        denovo_indel.append( variant_id)
                        w.writerow(newline)   
                    elif variant_pass and check_damage(variant, varclass):
                        w2.writerow(newline)              

        denovo_write.close()
        denovo2_write.close()
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


