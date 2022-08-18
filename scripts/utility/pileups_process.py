import argparse
import pandas as pd
import numpy as np
import os.path

'''Calculates pileup correlation and merges pileups for facets input'''


def rename_cols(col_name: str) -> str:
    '''
    Utility to rename columns in merging pileups
    :param: col_name: str, column name
    :returns: col_name: str, transformed name
    '''
    if '_x' in col_name:
        col_name = col_name[:-2]
    if '_y' in col_name:
        col_name = col_name[:-2].replace('1', '2')

    return col_name


# command line args
description_info = '''Calculates pileup correlation and merges pileups for facets input (pileups are not calculated together to be able to compare any sample if neeeded)'''

parser = argparse.ArgumentParser(
    description=(description_info),
    prog='pileups_process.py',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('patient_folder', type=str,
                    help='Patient folder to read&write csv')
args = parser.parse_args()
print("PROGRAM ARGUMENTS: ", vars(args))

patient_folder = args.patient_folder

# read pileups

pileups = []
for i in ['tumor-pileup.gz', 'normal-pileup.gz']:
    pileup = pd.read_csv(os.path.join(patient_folder, i))
    pileup = pileup[pileup['Alt'] != '.']

    pileups.append(pileup)

# merge for facets input

merged = pileups[1].merge(pileups[0], on=['Chromosome', 'Position', 'Ref', 'Alt'])

merged = merged.rename(mapper=rename_cols, axis=1).dropna()
merged.to_csv(os.path.join(patient_folder, 'merged-pileup.gz'), index=None)

# calculatiing correlation part

for idx, pileup in enumerate(pileups):

    pileup.index = pileup['Chromosome'] + '_' + pileup['Position'].astype(str)
    # drop low depth positions
    pileup = pileup[(pileup['File1R'] + pileup['File1A']) > 50]
    pileups[idx] = pileup.assign(VAF=pileup['File1A']/(pileup['File1R'] + pileup['File1A'] + 1))

# intersect on positions
cross_index = pileups[1].index.intersection(pileups[0].index)
pileups[0] = pileups[0].loc[cross_index]
pileups[1] = pileups[1].loc[cross_index]

# correlation

corr = np.corrcoef(pileups[0]['VAF'], pileups[1]['VAF'])[0][1]

with open(os.path.join(patient_folder, 'tumor-normal-correlation.txt'), 'w') as f:
    f.write(str(corr))
  