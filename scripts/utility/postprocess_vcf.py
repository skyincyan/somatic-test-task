import argparse
import pandas as pd
import gzip
import os.path
from pybedtools import BedTool

'''Adds custom annotations, custom filters to vcf. Filters by higher QC rates, tumor and normal depths, vafs, adds FFPE filter, adds total and minor copynumbers for each mutation'''

# filtering constants 
T_DEPTH_MIN = 10
N_DEPTH_MIN = 10
T_ALT_COUNTS_MIN = 10
N_VAF_MAX = 0.05
T_VAF_MIN = 0.05
FFPE_PARAM = 1
QSS_QSI_MIN = 20
EVS_MIN = 8

# used mappings
pass_dict = {True: 'PASS', False: 'LowQC'}  # filter mapping
nucleotide_dict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}  # count columns for vcf

# command line args
description_info = '''Adds custom annotations, custom filters to vcf. Filters by higher QC rates, tumor and normal depths, vafs, adds FFPE filter, adds total and minor copynumbers for each mutation'''

parser = argparse.ArgumentParser(
    description=(description_info),
    prog='postprocess_vcf.py',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('patient_folder', type=str, help='Patient folder to read&write csv')
parser.add_argument('preparation', type=str, help='How sample was prepared, FF or FFPE')

args = parser.parse_args()
print("PROGRAM ARGUMENTS: ", vars(args))

patient_folder = args.patient_folder
preparation = args.preparation

# read_vcf
vcf = pd.read_csv(os.path.join(patient_folder, 'strelka_somatic_merged_filtered_annotated.vcf.gz'), sep='\t', comment='#', header=None,
                     names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR'])

# read comments
comments = []
with gzip.open(os.path.join(patient_folder, 'strelka_somatic_merged_filtered_annotated.vcf.gz'), 'r') as f:
    a = f.readline()
    while (a[1:].decode().startswith('#')):
        comments.append(a)
        a = f.readline()
comments.append(a)

# extract info for filtering in each row and calculate needed stats

new_cols = []
for idx, row in vcf.iterrows():

    if row['INFO'].split(';')[1].startswith('QSS'):
        # snp counts
        refCounts = int(row['TUMOR'].split(':')[nucleotide_dict[row['REF']]].split(',')[0])
        altCounts = int(row['TUMOR'].split(':')[nucleotide_dict[row['ALT']]].split(',')[0])

        nrefCounts = int(row['NORMAL'].split(':')[nucleotide_dict[row['REF']]].split(',')[0])
        naltCounts = int(row['NORMAL'].split(':')[nucleotide_dict[row['ALT']]].split(',')[0])
    else:
        # indel counts
        refCounts = int(row['TUMOR'].split(':')[2].split(',')[0])
        altCounts = int(row['TUMOR'].split(':')[3].split(',')[0])

        nrefCounts = int(row['NORMAL'].split(':')[2].split(',')[0])
        naltCounts = int(row['NORMAL'].split(':')[3].split(',')[0])

    t_depth = int(row['TUMOR'].split(':')[0])
    n_depth = int(row['NORMAL'].split(':')[0])

    qs = int(row['INFO'].split(';')[1].split('=')[1])
    evs = float(row['INFO'].split(';')[11].split('=')[1])

    # calc vaf
    nvaf = naltCounts/(naltCounts+nrefCounts+1)
    tvaf = altCounts/(altCounts+refCounts+1)

    new_cols.append([altCounts, tvaf, t_depth, n_depth, nvaf, qs, evs])

vafs = pd.DataFrame(new_cols, columns=['t_alt_counts', 'Tumor_VAF', 't_depth', 'n_depth', 'Normal_VAF', 'QS', 'EVS'])

# base filters 

basic_thresholds = (vafs['t_alt_counts'] > T_ALT_COUNTS_MIN) & (vafs['t_depth'] > T_DEPTH_MIN) & (vafs['n_depth'] > N_DEPTH_MIN) & (vafs['Normal_VAF'] < N_VAF_MAX) & (vafs['Tumor_VAF'] > T_VAF_MIN)

higher_filter = (vafs['QS'] > QSS_QSI_MIN) & (vafs['EVS'] > EVS_MIN)

# add FFPE filter is needed and generate string for INFO 

if preparation == 'FFPE':
    ffpe_filter = (vafs.Tumor_VAF * vafs.t_alt_counts >= FFPE_PARAM)
    combined_filter = ffpe_filter & higher_filter & basic_thresholds
    filter_string = ';AFA=' + higher_filter.map(pass_dict) + '|' + ffpe_filter.map(pass_dict) + '|' + basic_thresholds.map(pass_dict) + '|' + combined_filter.map(pass_dict)
else:
    combined_filter = higher_filter & basic_thresholds
    filter_string = ';AFA=' + higher_filter.map(pass_dict) + '|PASS|' + basic_thresholds.map(pass_dict) + '|' + combined_filter.map(pass_dict)

# get segments to annotate mutations and prepare for bedtools
segments = pd.read_csv(os.path.join(patient_folder, 'segments.txt'), sep='\t')
segments_short = segments[['chrom', 'start', 'end', 'tcn.em', 'lcn.em']]

vcf_for_bed = pd.DataFrame({
    'chrom': vcf['CHROM'].map(lambda x: x[3:].replace('X', ' 23').replace('Y', '24')).astype(int), 
    'start': vcf['POS'], 'end': vcf['POS']+vcf['ALT'].map(len),
    'name': vcf.index
})   

# intersecting with bedtools
bed_file = BedTool.from_dataframe(segments_short)
intersection = BedTool.from_dataframe(vcf_for_bed).intersect(bed_file, wo=True).to_dataframe(header=None)

intersection.index = intersection['name']
intersection = intersection.reindex(vcf.index)

# add segments to INFO string
filter_string = filter_string + '|' + intersection['thickEnd'].astype(str) + '|' + intersection['itemRgb'].astype(str) 

vcf = vcf.assign(INFO=vcf['INFO']+filter_string)

# write result
with gzip.open(os.path.join(patient_folder, 'strelka_somatic_merged_filtered_postprocessed.vcf.gz'), 'w') as f:
    for c in comments[:-1]:
        f.write(c)
    f.write('##INFO=<ID=AFA,Number=.,Type=String,Description="Additional custom filtering and annotation. Format: HIGHER_CUTOFF|FFPE_FILTER|BASIC_FILTER|COMBINED_FILTER|Total_CN|Minor_CN">\n'.encode())
    f.write(comments[-1])

vcf.to_csv(os.path.join(patient_folder, 'strelka_somatic_merged_filtered_postprocessed.vcf.gz'), sep='\t', header=None, index=None, mode='a')
