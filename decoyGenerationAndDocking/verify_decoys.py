'''
This script cross references the DeepCoy generated decoys
to the DeepCoy selected decoys to check for missing decoys.
When found, it collates a new smi file containing the
original DeepCoy generated decoys for re-evaluation in
in order to produce adequate decoys for all actives in the
dataset.
'''

import os
from tqdm import tqdm
import pandas as pd
import shutil
import sys

def load_smiles_to_df(smiles_folder): # parses active and paired decoy smile into dataframe

    # make absolute paths for each smile file
    smiles_files = [f'{smiles_folder}{file}' for file in os.listdir(smiles_folder)]

    # set up lists for populating
    file_actives = list()
    file_decoys = list()

    # extract active and decoy from line of smile file
    print(f'Loading smiles from {smiles_folder}...')
    with tqdm(total=len(smiles_files)) as pbar:
        for file in smiles_files:
            with open(file, 'r') as text:
                lines = text.read()
                for line in lines.split('\n'):
                    try:
                        active, decoy = line.split(' ')[0], line.split(' ')[1]
                        file_actives.append(active)
                        file_decoys.append(decoy)
                    except:
                        pass
            pbar.update(1)

    # load actives and decoy smiles into df
    df = pd.DataFrame({'active':file_actives, 'decoys':file_decoys})

    return df

def output_smiles_for_re_evaluation(eval_output_folder, raw_output_folder): # find actives with less than 15 decoys and raw decoys to output file

    # load smiles
    eval_df = load_smiles_to_df(eval_output_folder)
    raw_df = load_smiles_to_df(raw_output_folder)

    # find actives with less than 15 selected decoys
    print('Checking for actives with less than 15 decoys...')
    re_evaluate = list()

    with tqdm(total=len(raw_df.active.unique())) as pbar:
        for active in raw_df.active.unique():
            post_eval_count = len(eval_df.loc[eval_df.active == active])
            if post_eval_count < 15:
                re_evaluate.append(active)
            pbar.update(1)

    # fetch generated decoys for actives missing evaluated decoys
    output_df = raw_df.loc[raw_df.active.isin(re_evaluate)]

    print(f'Found {len(output_df)/1000} actives with less than 15 decoys')

    print('Writing outfile...')
    # write to new file for re-evaluation by DeepCoy
    with open(f'/evaluation/re_evaluate_decoys.smi','a+') as final:
        for active, decoy in zip(list(output_df['active']), list(output_df['decoys'])):
            final.write(f'{active} {decoy}\n')
        final.close()

def main(): # run script using CLI

    output_smiles_for_re_evaluation('/final_decoys/','/screened_output/')

if __name__ == '__main__':
    main()
