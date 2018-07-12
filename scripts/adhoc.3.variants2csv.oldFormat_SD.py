#takes rare inherited vcf directory as input output summary csv
# python /home/local/ARCS/hq2130/Exome_Seq/scripts/python_script/adhoc.4.variants_to_csv.py -v vcf_inherited -o CHD_MedExomeKit_inherited0318.csv
import csv
import gzip
import os
from optparse import OptionParser
from collections import defaultdict

'''
takes denovo/inherited trio vcf folder as input, aggregate all variants into one csv without filtering
one variant per line


Usage example:   

    python $EXOMPY/adhoc.3.variants2csv.py -v $VCF_folder'/vcf_rare_denovo'  -o $denovo_csv
'''


def info_list(INFOstring, head):
    
    INFOcolumnList=INFOstring.split(";")
    INFOdict={}
    for element in INFOcolumnList:
        if '=' in element:
            FieldName,FieldValue=element.split('=',1)
            INFOdict[FieldName]=FieldValue           
    result = []
    for e in head:
        result.append(INFOdict.get(e,'NA'))
        
    return result


parser = OptionParser()
parser.add_option("-v", "--vcf", dest="VCFfolder",help="input VCF folder", metavar="VCFfolder")
parser.add_option('-o', '--output', dest='output',help='output csv name', metavar='output')
(options, args) = parser.parse_args()
vcf_folder = options.VCFfolder
output = options.output


print '-' * 50
print 'aggregate variants from vcf to csv'


head = set()
mutation = {}
mutation_counts = defaultdict(lambda: [0, 0])

for e in os.listdir(vcf_folder):
    if e.endswith('.vcf'):
        with open(vcf_folder +'/'+ e) as f1:
            proband, father_ID, mother_ID= e.rsplit('.vcf')[0].split('_')
            
            for line in f1:
                if line[0] =='#':
                    if line.startswith('##INFO=<ID='):
                        head.add(line.split('=')[2].split(',')[0])
                else:
                   
                    info = line.split()
                    chrom = info[0]
                    pos = info[1]
                    INFOstring = info[7]
                    format = info[8]

                    variant_id = '_'.join([proband, chrom, pos])
                    variant_id_2 = '_'.join([chrom, pos])

                    proband_geno = info[9].split(':')[0]
                    proband_info = info[9]
                    Father_geno = info[10].split(':')[0]
                    Father_info = info[10]
                    Mother_geno = info[11].split(':')[0]
                    Mother_info = info[11]
                    
                    parents_list = father_ID +'('+Father_info+')'+','+mother_ID +'('+ Mother_info +')'
                    mutation_counts[variant_id_2][0] += 1

                    if variant_id in mutation:
                        print 'same varaint in same proband again:' , variant_id
                    else:
                        mutation[variant_id] = [format, proband+'('+ proband_info +')', parents_list] + info[:8]
head = list(head)
head.append("Calt")
head.sort()
#print e+"\t"+output
with open(output,'wb') as f:
    w = csv.writer(f)
    w.writerow(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','FORMAT','proband','parents', 'het', 'homo'] + head)
    for variant_id, value in mutation.items():
        variant_id_2 = '_'.join(variant_id.split('_')[-2:])
        #print variant_id, value
        content =  value[3:-1] + value[:3] + mutation_counts[variant_id_2] + info_list(value[-1], head)
        w.writerow(content)

    
print '-' * 50
print 'aggregate variants to csv finished'
print "output to "+output
