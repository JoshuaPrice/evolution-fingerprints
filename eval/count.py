'''
Counts number of SNPs in each chromosome
Prints the list of chromosomes sorted by SNP count
'''

import sys
import csv

import numpy

assert len(sys.argv) > 1, 'Usage: count.py list of csv files to count'
files = sys.argv[1:]

chrom_count = {}
for filename in files:
    with open(filename, newline='') as rsfile:
        reader = csv.reader(rsfile, delimiter=',')
        for row in reader:
            chrom = int(row[0])
            if chrom in chrom_count.keys():
                chrom_count[chrom] += 1
            else:
                chrom_count[chrom] = 1

for key, value in sorted(chrom_count.items(), key=lambda x: x[0]):
    print(key, ":", value)
