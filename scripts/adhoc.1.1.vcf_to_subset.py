#  takes a vcf file, filtered rare exon variants and output to _rare.vcf

import os
from optparse import OptionParser


usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-s", "--subset", dest="SUBfile",help="input samples file", metavar="SUBfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")

(options, args) = parser.parse_args()
vcf_file = os.path.abspath(options.VCFfile)
subset_file = os.path.abspath(options.SUBfile)
out_file = options.OutputFileName

with open(vcf_file, 'r') as f:
    for line in f: 
        if line[:6]=='#CHROM': 
            samples = line.strip().split('\t')
            break

case_index=[]
with open(subset_file) as f:
    for line in f:
        case = line.strip().split()[0]
        if case in samples:
            case_index.append(samples.index(case))

print len(case_index)

ChrCount=0
with open(vcf_file, 'r') as f:
    fw = open(out_file, 'w')
    for line in f: 
        if line[:2]=='##': # get and write head
            fw.write(line)
        else: 
            data = line.strip().split('\t')
            # verbose  
            ChrPresent = data[0]
            if ChrCount != ChrPresent:
                print "Chromosome "+str(ChrPresent)
                ChrCount=ChrPresent

            subset = []
            for index in case_index:
                subset.append(data[index])
            if '0/0' in subset[0] or  './.' in subset[0]:
                continue
            
            fw.write('\t'.join(data[:9] + subset) + '\n')
           
    fw.close()

            
