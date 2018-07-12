#!/usr/bin/python
from __future__ import division
import csv
import gzip
import shutil
import subprocess
import os
import sys
import math
from optparse import OptionParser
from collections import defaultdict



#vcfname = 'gvcf.rawvariants.hardfiltered.vcf'
#newvcf = 'PCGC_mappability.vcf'
vcfname = sys.argv[1]
newvcf = sys.argv[2]
print vcfname
print newvcf
with open(vcfname) as f:
    fw = open('temp.txt','w')
    for line in f:
        if not line.startswith('#'):
            lst = line.split()
            fw.write('\t'.join(['chr'+lst[0],lst[1],lst[1],'\n']))
    fw.close()
 
print 'bedtool rmsk'  
dir_map = '/home/local/ARCS/hq2130/Exome_Seq/resources/mappability/'
cmd = 'bedtools intersect -a '+dir_map+'hg19.rmsk.bed -b temp.txt'
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
repeat = p.communicate()[0]

rep = {}
for hit in repeat.split('\n'):
    if len(  hit.split()) > 0:
        chrom, start, end, repname = hit.split()[:4]
        rep['_'.join([chrom, start, end])] = repname

print 'bedtool mappability'  
cmd = 'bedtools intersect -a '+dir_map+'wgEncodeCrgMapabilityAlign75mer.bedGraph -b temp.txt'
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
repeat = p.communicate()[0]

map_75bp = {}
for hit in repeat.split('\n'):
    if len(  hit.split()) > 0:
        chrom, start, end, mappability = hit.split()[:4]
        map_75bp['_'.join([chrom, start, end])] = mappability

print 'write to new vcf' 
with open(vcfname) as f:
    fw = open(newvcf,'w')
    for line in f:
        if not line.startswith('#'):
            lst = line.split()
            chrom, start, end = 'chr' + lst[0], lst[1], lst[1]
            variant_id = '_'.join([chrom, start, end])
            
            region,  mappability = 'No', 'NA'
            if variant_id in rep:
                region = 'Yes'
            if variant_id in map_75bp:
                mappability = map_75bp[variant_id]
            lst[7] = lst[7] + ';' + ';'.join(['RepeatMask='+region, '75bpMap='+mappability])
            
            fw.write('\t'.join(lst)+'\n')
        else:
            fw.write(line)
    fw.close()
    
    

os.remove('temp.txt')


