'''
This script checks how many decoys are present for each active
in the user supplied active pdbqts folder, produces a summary
csv of decoy counts for each active, and offers to make a copy
of the dataset with only actives that have more than a user
defined threshold number of decoys. Also offers to copy the
decoys for these actives only into te same folder.
'''

import os
from tqdm import tqdm
import pandas as pd
import shutil
import sys

def remove_empty_folders(decoy_pdbqts): # delete empty folders in docked decoys folder due to GWOVina errors

    print('Removing empty folders...')
    decoy_folders = [f'{decoy_pdbqts}{folder}' for folder in os.listdir(decoy_pdbqts)]
    with tqdm(total=len(decoy_folders)) as pbar:
        for folder in decoy_folders:
            if len(os.listdir(folder)) != 2:
                os.rmdir(folder)
            pbar.update(1)

def count_decoys_for_structures(decoy_pdbqts, active_pdbqts): # check how many decoys have been generated for each active

    # set up for populating
    decoy_active_counts = dict()

    # make a directory for problem decoys
    parent_path = decoy_pdbqts.split('/')[:len(decoy_pdbqts.split('/'))-2]

    parent_path = os.path.join(*parent_path)

    if not os.path.isdir(f'/{parent_path}/problem_decoys'):
        os.mkdir(f'/{parent_path}/problem_decoys')

    # loop through structures and count how many decoy files contain pdb code string
    print('Counting Decoys...')
    with tqdm(total=len(os.listdir(active_pdbqts))) as pbar:
        for pdb_code in os.listdir(active_pdbqts):

            # find decoy folders for that structure
            decoys_with_pdb_code = [folder for folder in os.listdir(decoy_pdbqts) if pdb_code in folder]

            # count successfully docked decoys
            num_decoys = len(decoys_with_pdb_code)
            decoy_active_counts[pdb_code] = num_decoys
            pbar.update(1)

    # turn counts into a df
    decoy_summary_df = pd.DataFrame(decoy_active_counts.items(), columns=['Active_PDBCode','Num_Decoys'])
    decoy_summary_df.to_csv(os.path.join(os.path.dirname(__file__), 'actives_decoy_summary.csv'), index=False)

    return decoy_summary_df, parent_path

def summarise_missing(decoy_thresh, decoy_summary_df, active_pdbqts, decoy_pdbqts, parent_path): # print summary of decoy counts and make new dataset excluding actives with no decoys

    # summary
    print('Decoy summary:')
    print(decoy_summary_df['Num_Decoys'].value_counts())
    print(f'Found {len(decoy_summary_df.loc[decoy_summary_df.Num_Decoys < decoy_thresh])} structures with less than {decoy_thresh} decoys present:')
    no_decoys = list(decoy_summary_df.loc[decoy_summary_df.Num_Decoys < decoy_thresh]['Active_PDBCode'])
    print(f'{no_decoys}')

    # option to make new dataset with only structures that have decoys
    answer = input(f'Make copy of dataset without these structures under {active_pdbqts[:len(active_pdbqts) - 1]}_with_decoys? (y/n)')
    if answer == 'y':
        structures_with_decoys = [structure for structure in os.listdir(active_pdbqts) if structure not in no_decoys]
        os.mkdir(f'{active_pdbqts[:len(active_pdbqts) - 1]}_with_decoys')
        print('Copying structures...')
        with tqdm(total=len(structures_with_decoys)) as pbar:
            for structure in structures_with_decoys:
                shutil.copytree(f'{active_pdbqts}{structure}', f'{active_pdbqts[:len(active_pdbqts) - 1]}_with_decoys/{structure}')
                pbar.update(1)

        # option to copy decoys for these structures to this dataset also
        include = input('Copy decoys into this dataset? (y/n)')
        if include == 'y':
            print('Copying decoys...')
            with tqdm(total=len(structures_with_decoys)) as pbar:
                for structure in structures_with_decoys:

                    # get folder names of all decoys for that active
                    decoys_to_copy = [decoy for decoy in os.listdir(decoy_pdbqts) if decoy[:4] in structure]

                    # temporarily change naming convention to from decoy_1 to decoy_01 for sorting to ensure best 10 decoys are selected
                    decoys_to_copy = sorted([f'{decoy[:11]}0{decoy[11:]}' if len(decoy[11:]) == 1 else decoy for decoy in decoys_to_copy])[:decoy_thresh]
                    for decoy in decoys_to_copy:
                        
                        # amend naming convention in filepath for copying
                        if decoy[11] == '0':
                            decoy = f'{decoy[:11]}{decoy[12]}'
                            shutil.copytree(f'{decoy_pdbqts}{decoy}', f'{active_pdbqts[:len(active_pdbqts) - 1]}_with_decoys/{decoy}')
                        else:
                            shutil.copytree(f'{decoy_pdbqts}{decoy}', f'{active_pdbqts[:len(active_pdbqts) - 1]}_with_decoys/{decoy}')
                    pbar.update(1)
        else:
            pass

def parse_args(args): # parse CLI user inputs

    decoy_pdbqts = args[args.index('-decoys') + 1]

    active_pdbqts = args[args.index('-pdbqts') + 1]

    decoy_thresh = int(args[args.index('-decoy_thresh') + 1])

    return decoy_pdbqts, active_pdbqts, decoy_thresh

def main(): # run script using CLI

    decoy_pdbqts, active_pdbqts, decoy_thresh = parse_args(sys.argv)

    # remove empty decoy folders
    remove_empty_folders(decoy_pdbqts)

    # count decoys into df
    decoy_summary_df, parent_path = count_decoys_for_structures(decoy_pdbqts, active_pdbqts)

    # summarise and make new dataset
    summarise_missing(decoy_thresh, decoy_summary_df, active_pdbqts, decoy_pdbqts, parent_path)

if __name__ == '__main__':
    main()
