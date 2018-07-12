

import gzip
import os
from optparse import OptionParser
from copy import deepcopy

'''
takes vcf and ped as input, seperate the vcf by trios(proband and parents are in the vcf), 
exclude varaints are all missing in the trio
example:    python adhoc.1.VcfTotrio.py -v PCGC_target3_0223.vcf -p CHD_MedExomeKit.ped

to do:
trio_GT not in {['0/0','0/0','0/0'],  ['./.','./.','./.']

add logging part 

'''

def info_split(INFOstring):
    INFOcolumnList=INFOstring.split(";")
    INFOdict={}
   
    for element in INFOcolumnList:
        if '=' in element:
            FieldName,FieldValue=element.split('=',1)
            INFOdict[FieldName]=FieldValue
    return INFOdict

def check_rare_exon(INFOdict,index):# find variant with propability < 0.1%
    threshold = 1 / 100
    CodingCodes={'splicing', 'exonic', 'exonic,splicing'}
    # Get values from INFO fields
   # if multiallele : 
    Astr="change_INFO=1"
    for key, value in INFOdict.iteritems():
        tempstr=value.split(',')
       # print key+'\t'+value;
        if len(tempstr) >1  and index <= len(tempstr):
            INFOdict[key]=tempstr[index-1];
        if index > len(tempstr) and len(tempstr) >1 :
           print value +tempstr[0];
        Astr=Astr+";"+key+"="+INFOdict[key]
    #return Astr
    INFOstring=Astr  
#    laf= float(INFOdict.get("AC"))/float(INFOdict.get("AN"))    
#    KGscore= max(map(float,[ rate  for rate in INFOdict.get("1KGfreq",str(threshold)).split(',') if rate != '.' ])+[0])
#    ESPscore= max(map(float,[ rate  for rate in INFOdict.get('ESPfreq',str(threshold)).split(',') if rate != '.' ])+[0])
#    ExAcscore= max(map(float,[ rate  for rate in INFOdict.get('ExACfreq',str(threshold)).split(',') if rate != '.' ])+[0])
    MutationFunct=str(INFOdict.get('VarFunc','None').split(',')[0])
    #print str(KGscore)+"\t"+str(ESPscore)+"\t"+str(ExAcscore)+'\t'+MutationFunct
    if MutationFunct == 'None':
        print 'strange, No variant function result'
        return ""
    if MutationFunct in CodingCodes :
         return INFOstring
#    if  KGscore <= threshold and ESPscore <= threshold and  ExAcscore <= threshold and MutationFunct in CodingCodes and laf<0.02 :
#        return INFOstring
#    else : 
#        return ""




# Basic Input Files, required   

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-p", "--ped", dest="PEDfile",help="input PED file", metavar="PEDfile")

(options, args) = parser.parse_args()
vcf_name = os.path.abspath(options.VCFfile)
ped_name = os.path.abspath(options.PEDfile)

# make vcf trio dir
trio_dir = os.path.dirname(ped_name)  
#trio_dir = os.getcwd()  
print trio_dir+" "+ped_name

trio_dir_name = trio_dir+'/vcf_trio/'
if not os.path.exists(trio_dir_name):
    os.makedirs(trio_dir_name)
    
if vcf_name.endswith(".gz"):
    f=gzip.open(vcf_name, 'r')
    #with gzip.open(vcf_name, 'r') as f:
else:
    f=open(vcf_name, 'r')
    #with open(vcf_name, 'r') as f:
for line in f: 
    if line[:6]=='#CHROM': 
        samples = line.split()
        break

# ulimit -n 2048 
sample_index={}
included_sample = set()
name_exceptions=["X","-",".","x"]
with open(ped_name) as f:
    for line in f:
        if not line.startswith("#") :
            family, proband,father,mother = line.strip().split('\t')[:4]
            # print family, proband,father,mother
        
            if proband in samples and father in samples and mother in samples:
                sample_index[proband]=[father,mother,samples.index(proband),samples.index(father),samples.index(mother)]
                included_sample.add(father)
                included_sample.add(mother)
                included_sample.add(proband)
            else:
                if proband in samples:
                    if father not in samples and mother not in samples :
                         sample_index[proband]=[father,mother,samples.index(proband),-1,-1]
                         included_sample.add(proband)
                    else:
                         if father not in samples:
                             sample_index[proband]=[father,mother,samples.index(proband),-1,samples.index(mother)]
                             #included_sample.add(father)
                             included_sample.add(mother)
                             included_sample.add(proband)
                         else:
                             sample_index[proband]=[father,mother,samples.index(proband),samples.index(father),-1]
                             included_sample.add(father)
                        #included_sample.add(mother)
                             included_sample.add(proband)
print '-' * 50

#print 'non-trio samples:', set(samples) - included_sample
print 'trios to seperate:',len(sample_index)


print 'seperate vcf started'
ChrCount=0
if vcf_name.endswith(".gz"):
    f=gzip.open(vcf_name, 'r')
else:
    f=open(vcf_name, 'r')
#with open(vcf_name, 'r') as f:
head = []
for line in f: 
    if line[0]=='#': # get and write head
        if line[1]=='#':
            # format and info line
            head.append(line)
          #  print(line)
        else: #
            all_samples = line.strip().split('\t')
            head_f = all_samples[:9] # CHROM ... FORMAT
 
            f = [open(trio_dir_name+ '__'.join([proband, sample_index[proband][0], sample_index[proband][1] ]) + '.vcf', "w") for proband in sample_index]
            n = len(f)
            for i in range(n):
                proband = f[i].name.split('/')[-1].split('__')[0]
                proband_index, father_index, mother_index = sample_index[proband][-3:]
               # print proband+"\t"+f[i].name+"\t"+ str(proband_index)
                fstr="-"
                mstr="-"
                if father_index != -1 : 
                    fstr=all_samples[father_index]
                if mother_index != -1 :
                    mstr=all_samples[mother_index]    
                        
                trio = [all_samples[proband_index],fstr,mstr]
                f[i].write(''.join(head))
                f[i].write('\t'.join(head_f + trio)+'\n')   
    else: 
        data = line.strip().split('\t')

        # verbose  
        ChrPresent = data[0]
        if ChrCount != ChrPresent:
            print "Chromosome "+str(ChrPresent)
            ChrCount=ChrPresent
       
        INFOdict=info_split(data[7])
        if int(INFOdict['AN']) < 2*(len(data)-9)*0.9 : continue
        alts=data[4].split(',')
        for i in range(n):
            proband = f[i].name.split('/')[-1].split('__')[0]
         #   print proband+"\t"+f[i].name+"\t"+ str(proband_index)
            proband_index, father_index, mother_index = sample_index[proband][-3:]
        #    print proband+"\t"+f[i].name+"\t"+ str(proband_index)
            fgt="./."
            mgt="./."
            if father_index != -1 :
                fgt=data[father_index]
            if mother_index != -1 :
                mgt=data[mother_index]
                
            trio = [data[proband_index],fgt,mgt]
            infostr=data[7]
            alt=data[4]
            # we replace ./. to 0/0 so that we can exclude sth like [0/0, ./., 0/0 ], not used anymore
            pgt0=data[proband_index].split(':')[0].replace('.','0')
           # print str(proband_index)+"\t"+pgt0
            cor=max([int(numeric_string) for numeric_string in pgt0.split('/')])
           # print pgt0+"\t"+str(cor)+"\t"+str(cor is 0)+"\t"+str(len(alts))
            if father_index != -1 :
                fgt=data[father_index].split(':')[0].replace('.','0')
                if cor < max([int(numeric_string) for numeric_string in fgt.split('/')]) :
                    cor=0;
            if mother_index != -1 :
                mgt=data[mother_index].split(':')[0].replace('.','0')
                if cor < max([int(numeric_string) for numeric_string in mgt.split('/')])  : 
                    cor=0;
            trio_GT = [pgt0,fgt,mgt]
          
            if cor > 0 and len(alts) >1 :
                alt=alts[cor-1]
          # print '\t'.join(data[:3])+"\t"+alt+'\t'
            
           # print trio_GT[0] +"\t"+ str (cor) 
            if trio_GT != ['./.','./.','./.'] and trio_GT[0] !="0/0" and cor>0  and len(alts)<4 :
                INFOdict2=deepcopy(INFOdict)
                infostr=check_rare_exon(INFOdict2,cor);
                if infostr :
                    wstr='\t'.join(data[:4])+'\t'+alt+'\t'+'\t'.join(data[5:7])+'\t'+"Calt="+str(len(alts))+";"+infostr+'\t'+'\t'.join(data[8:9])+'\t'+'\t'.join(trio)+'\n'
                    #    print str(cor)+'\t'+data[3]+"\t"+data[4]+"\t"+alt
                 #   print str(cor)+'\t'+trio_GT[0]+'\t'+data[4]
                  #  print(wstr)
                    f[i].write(wstr)

for fh in f:
    fh.close()

print 'seperate vcf finished'
print '-' * 50
print ''
