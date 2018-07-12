#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import string
import tempfile
import argparse
import subprocess


def lift_over(chrom, pos, strand):
    lift_over_exec ="/home/local/ARCS/nz2274/Application/liftover/liftOver"
    map_chain ="/home/local/ARCS/nz2274/Application/liftover/hg19ToHg38.over.nochr.chain.gz" #os.path.join(os.path.dirname(__file__), 'chain', chain + '.over.chain.gz')

    with tempfile.NamedTemporaryFile() as query:
        q = '\t'.join([ chrom,  # chrom in .chain startswith `chr` #'chr' + chrom if not chrom.startswith('chr') else 
                       str(pos - 1),  # 1-based to 0-based
                       str(pos),      # 1-based to 0-based
                       '.',           # ignored
                       '0',           # ignored
                       strand])
        query.write(q)
        query.write('\n')
        query.seek(0)

        with tempfile.NamedTemporaryFile() as mapped, tempfile.NamedTemporaryFile() as unmapped:
            # $ ./liftOver <oldFile> <map.chain> <newFile> <unMapped>
            cmd = [lift_over_exec, query.name, map_chain, mapped.name, unmapped.name]
            # print >>sys.stderr, '[DEBUG]', cmd
            subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]

            mapped.seek(0)
            for line in mapped:
                record = line.strip().split('\t')
                return (record[0].replace('chr', ''), record[2], record[5])

def reverse_complement(x):
    reverse_complement_map = string.maketrans('ATGC', 'TACG')
    return x.translate(reverse_complement_map)

def _main():
    parser = argparse.ArgumentParser()
    #parser.add_argument('--chain', choices=['hg18ToHg19', 'hg18ToHg38', 'hg19ToHg38'], required=True)
    parser.add_argument('--format', choices=['vcf'], default='vcf')
    args = parser.parse_args()


    for line in sys.stdin:
        if args.format == 'vcf':
            if line.startswith('#'):
                print line,
            else:
                record = line.strip().split('\t')
                new = lift_over(record[0],
                                int(record[1]),  # 1-based pos
                                '+'             # vcf is always forward to reference
                                )
                if new:
                    new_chrom, new_pos, new_strand = new
                    ref = record[3] if new_strand == '+' else reverse_complement(record[3])
                    alt = record[4] if new_strand == '+' else reverse_complement(record[4])
                    print '\t'.join([new_chrom, new_pos, record[2], ref, alt] + record[5:])


if __name__ == '__main__':
    _main()
