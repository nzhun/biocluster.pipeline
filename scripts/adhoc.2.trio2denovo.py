# takes vcf directory, output new vcf with only denovo mutations
# general results, need to refine MAF, ref/(ref+alt)
# rare, pass, denovo
import os
import re
from optparse import OptionParser

'''
takes trio vcf folder input, call rare exon de novo variants and rare exon inherited varaiants in the trio

example:   

To do:

add flag to do rara or exon
'''


def check_rare_exon(INFOstring,index):# find variant with propability < 0.1%
    threshold = 0.1 / 100
    CodingCodes={'splicing', 'exonic', 'exonic,splicing'}
    
    INFOcolumnList=INFOstring.split(";")
    INFOdict={}
    multiallele=0;
    for element in INFOcolumnList:
        if '=' in element:
            FieldName,FieldValue=element.split('=',1)    
            if FieldName !="VQSR" and FieldName in INFOdict :
                FieldValue=INFOdict.get(FieldName)+","+FieldValue
                multiallele=1
              #  print FieldValue
            INFOdict[FieldName]=FieldValue
    # Get values from INFO fields
    if multiallele : 
        Astr="change_INFO=1"
        for key, value in INFOdict.iteritems():
            tempstr=value.split(',')
            if len(tempstr) >1 :
              #  print tempstr[0]+"  --- "+tempstr[1]+"---"+tempstr[index-1] +"---"+str(index)
                INFOdict[key]=tempstr[index-1];
            
            Astr=Astr+";"+key+"="+INFOdict[key]
        INFOstring=Astr  
         
    KGscore= max(map(float,[ rate  for rate in INFOdict.get("1000g2015aug_all",str(threshold)).split(',') if rate != '.' ])+[0])
    ESPscore= max(map(float,[ rate  for rate in INFOdict.get('esp6500siv2_all',str(threshold)).split(',') if rate != '.' ])+[0])
    ExAcscore= max(map(float,[ rate  for rate in INFOdict.get('ExAC_ALL',str(threshold)).split(',') if rate != '.' ])+[0])
    MutationFunct=str(INFOdict.get('Func.refGene','None').split(',')[0])
    if MutationFunct == 'None':
        print 'strange, No variant function result'
    if  KGscore <= threshold and ESPscore <= threshold and  ExAcscore <= threshold and MutationFunct in CodingCodes :
        return INFOstring
    else : 
        return ""


def check_denovo(trio):

    proband= trio[0].split(':')
    Father = trio[1].split(':')
    Mother = trio[2].split(':')
    proband_GT = proband[0]
    Father_GT = Father[0]
    Mother_GT = Mother[0]
    
    if proband_GT in {'./.','0/0'}: #proband genotype missing or not denovo
        return False
    else:
        if Father_GT == '0/0' and Mother_GT == '0/0':
            return True
        else:
            return False

def get_GQ(dicts):
    GQ = dicts.get('GQ','0')
    if GQ.isdigit():
        return int(GQ)
    else:
        return 0

def check_inherited(trio, name):# trio =[s,f,m]

    proband = dict(zip(name.split(':'),trio[0].split(':')))
    Father = dict(zip(name.split(':'),trio[1].split(':')))
    Mother = dict(zip(name.split(':'),trio[2].split(':')))
    proband_GT = proband['GT']
    proband_GQ = get_GQ(proband)
    Father_GT = Father['GT']
    Father_GQ = get_GQ(Father)
    Mother_GT = Mother['GT']
    Mother_GQ = get_GQ(Mother)

    GQ_threshold = 30
    GQ_threshold2 = 60
    GQ_threshold3 = 80
    
    # not inherited or unknown
    inherited = False
    if proband_GT  in {'./.', '0/0'} : 
        inherited = False
    else:
        # homozygous
        if proband_GT[0] == proband_GT[-1]: 
            # find in both parents
            if proband_GT[-1] in Father_GT and proband_GT[-1] in Mother_GT:
                if proband_GQ >= GQ_threshold and Father_GQ >= GQ_threshold and Mother_GQ >= GQ_threshold:
                    inherited = True
            # found in father, unknown in mother
            elif proband_GT[-1] in Father_GT and Mother_GT in {'./.'}:
                if proband_GQ >= GQ_threshold2 and Father_GQ >= GQ_threshold:
                    inherited = True
            elif proband_GT[-1] in Mother_GT and Father_GT in {'./.'}:
                if proband_GQ >= GQ_threshold2 and Mother_GQ >= GQ_threshold:
                    inherited = True
            # both unknown
            else:
                if proband_GQ >= GQ_threshold3 and Mother_GT == './.' and  Father_GT == './.' :
                    inherited = True
        # het
        else:
            if proband_GT[-1] in Father_GT and proband_GT[-1] in Mother_GT:
                if proband_GQ >= GQ_threshold and max(Father_GQ, Mother_GQ) >= GQ_threshold:
                    inherited = True
            elif proband_GT[-1] in Father_GT:
                if proband_GQ >= GQ_threshold and Father_GQ >= GQ_threshold:
                    inherited = True
            elif proband_GT[-1] in Mother_GT:
                if proband_GQ >= GQ_threshold and Mother_GT >= GQ_threshold:
                    inherited = True
            else:
                # both unknown
                if proband_GQ >= GQ_threshold2 and ( Mother_GT == './.' or  Father_GT == './.' ):
                    inherited = True

    return inherited

parser = OptionParser()
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF folder", metavar="VCFfile")
(options, args) = parser.parse_args()
vcf_name = options.VCFfile

# make vcf trio dir
trio_dir = os.path.abspath(vcf_name)
denovo_dir = '/'.join(trio_dir.split('/')[:-1])+'/vcf_rare_denovo/'
if not os.path.exists(denovo_dir):
    os.makedirs(denovo_dir)
inherited_dir = '/'.join(trio_dir.split('/')[:-1])+'/vcf_rare_inherited/'
if not os.path.exists(inherited_dir):
    os.makedirs(inherited_dir)

print '-' * 50
print 'filter de novo/inherited variants from vcf files into filtered vcf files'

idx = 0 
for fname in os.listdir(trio_dir):
    if fname.endswith('.vcf'):
        idx += 1
        #print 'processed {0} trios, now process {1}'.format(idx, fname)
        f1 = open(trio_dir+'/'+fname)
        fw1 = open(denovo_dir+'/'+fname,'w')
        fw2 = open(inherited_dir+'/'+fname ,'w')
        for line in f1:
            if line[0]=='#': # get and write head
                fw1.write(line)
                fw2.write(line)
            else: #write rare-deleterious-inherited mutations           
                data =  line.split()
                alts=data[4]
                INFOstring = data[7]
                trio_gt= data[9];
                ref,alt=re.split(':|/',trio_gt)[0:2];
                if ref =="." or alt=="." :
                    continue
                index=max(int(ref), int(alt))
                if index==0 :
                    continue
                INFOstring=check_rare_exon(INFOstring,index)
                if INFOstring: #population freq later
                    data[7]=INFOstring
                    if len(alts)>1 and index <len(alts.split(',')) :
                        data[4]=alts.split(',')[index-1]
                    trio = data[9:12]       
                    if check_denovo(trio):
                        fw1.write("\t".join(data[:9]+trio)+'\n') 
                    if check_inherited(trio,data[8]):
                        fw2.write("\t".join(data[:9]+trio)+'\n')         
        fw1.close()
        fw2.close() 

print 'call de novo/inherited from vcf finished'
print '-' * 50
print ''
