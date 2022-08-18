import argparse
import pandas as pd
import numpy as np
import gzip
import os.path

'''Script to finish snp and indel merge and create filtered vcf for annotating'''

def merge_columns(item1, item2):
    '''
    Utility to merge data in tumor/normal columns
    :param: item1: str, value col1
    :param: item2: str, value col2
    :returns str, final value
    '''
    if item1.startswith('.'):
        return item2
    else:
        return item1

# command line args 
description_info = '''Finishes vcf merge and prepares filtered vcf by FILTER == 'PASS' to annotate. (not annotating lowEvs and lowDepth). Works with preset names in supplied patient folder'''

parser = argparse.ArgumentParser(
    description=(description_info),
    prog='clear_merge_n_filter.py',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('patient_folder', type=str, help='Patient folder to read&write csv')
args = parser.parse_args()
print("PROGRAM ARGUMENTS: ", vars(args))

patient_folder = args.patient_folder

# read_vcf
vcf = pd.read_csv(os.path.join(patient_folder, 'merged_vcf.vcf.gz'), sep='\t', comment='#', header=None,
                 names=['CHROM', 'POS', 'ID','REF', 'ALT', 'QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR', 'NORMAL2','TUMOR2'])

# prepare comments
comments = []
with gzip.open(os.path.join(patient_folder, 'merged_vcf.vcf.gz'), 'r') as f:
    a = f.readline()
    while (a.decode().startswith('#')):
        comments.append(a)
        a = f.readline()

# modify header
comments[-1] = comments[-1][:-18] + "\n".encode()

# merge cols 
vcf = vcf.assign(TUMOR=vcf.apply(lambda x: merge_columns(x['TUMOR'], x['TUMOR2']), axis=1), NORMAL=vcf.apply(lambda x: merge_columns(x['NORMAL'], x['NORMAL2']), axis=1))
vcf = vcf.drop(['NORMAL2', 'TUMOR2'], axis=1)

# write unfiltered
with gzip.open(os.path.join(patient_folder, 'strelka_somatic_merged_unfiltered.vcf.gz'), 'w') as f:
    for c in comments: 
        f.write(c)
vcf.to_csv(os.path.join(patient_folder, 'strelka_somatic_merged_unfiltered.vcf.gz'), sep='\t', header=None, index=None, mode='a')

# not going to annotate LowEvs, very low quality anyway
vcf = vcf[vcf['FILTER'] == 'PASS']

# write filtered
with gzip.open(os.path.join(patient_folder, 'strelka_somatic_merged_filtered.vcf.gz'), 'w') as f:
    for c in comments: 
        f.write(c)
vcf.to_csv(os.path.join(patient_folder, 'strelka_somatic_merged_filtered.vcf.gz'), sep='\t', header=None, index=None, mode='a')