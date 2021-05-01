'''
This script splits the data ensuring there is an equal
distribution of both structure resolution and good:bad
binder ratio. It also ensures that decoys and actives
for the same structure do not span datasets to avoid
look ahead bias
'''

from pypdb import *
import pprint
import pandas as pd
from tqdm import tqdm
from sklearn.model_selection import StratifiedShuffleSplit
import matplotlib.pyplot as plt
import concurrent.futures
import sys
from warnings import filterwarnings
filterwarnings("ignore")

# set up for populating
resolution_dict = dict()

def get_resolution(pdb): # build dictionary of resolutions from PDB API for list of pdb codes

    resolution = get_info(pdb)['pdbx_vrpt_summary']['pdbresolution']

    # add resolution to dictionary
    resolution_dict[pdb] = float(resolution)

def batch_fetch_resolution(pdb_codes): # wrapper for multiprocessing convert_to_pdbqt function
    threads = min(48, len(pdb_codes))

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(get_resolution, pdb_codes), total=len(pdb_codes)))

def parse_args(args): # parse CLI user inputs

    full_data = args[args.index('-data') + 1]

    output_dir = args[args.index('-out') + 1]

    test_frac = float(args[args.index('-test_size') + 1])

    return full_data, output_dir, test_frac

def main(): # run script using CLI

    full_data, output_dir, test_frac = parse_args(sys.argv)

    data = pd.read_hdf(full_data, key='df', mode='r') # read in the full dataset

    actives_only = data.loc[~data.PDBCode.str.contains('decoy')] # drop the decoys in order to split the dataset avoiding decoys from the same structure spanning train and test sets

    pdbs_in_data = list(actives_only['PDBCode']) # master list of pdb codes for fetching resolution

    # fetch resolution of active structures
    print('Fetching structure resolutions for actives...')
    batch_fetch_resolution(pdbs_in_data)

    actives_only['Resolution'] = actives_only['PDBCode'].map(resolution_dict) # add resolution as column in df

    # bin resolution and binding data values to perform stratified split
    actives_only['Binding_Data_in_uM_Binned'] = pd.qcut(actives_only['Binding_Data_in_uM'], q=6, labels=[1, 2, 3, 4, 5, 6])
    actives_only['Resolution_Binned'] = pd.qcut(actives_only['Resolution'], q=6, labels=[1, 2, 3, 4, 5, 6])

    # split data based on label, resolution and binding data frequency
    splitter = StratifiedShuffleSplit(n_splits = 1, test_size=test_frac, random_state=12)
    for train_index, test_index in splitter.split(actives_only, actives_only[['Label', 'Binding_Data_in_uM_Binned', 'Resolution_Binned']]):
        training_set = actives_only.iloc[train_index]
        test_set = actives_only.iloc[test_index]

    # split data based on returned fair indexes for splitting
    training_pdb_codes = list(training_set['PDBCode'])
    test_pdb_codes = list(test_set['PDBCode'])

    # add active specific decoys back into training and test sets
    full_training_codes = [code for code in list(data['PDBCode']) if code[:4] in training_pdb_codes]
    full_test_codes = [code for code in list(data['PDBCode']) if code[:4] in test_pdb_codes]
    full_training_df = data.loc[data['PDBCode'].isin(full_training_codes)]
    full_test_df = data.loc[data['PDBCode'].isin(full_test_codes)]

    # summarise split for user
    print(f'All Data Length: {len(data)}')
    print(f'Training Data Length: {len(full_training_df)}')
    print(f'Test Data Length: {len(full_test_df)}')
    print(f'Combined Length: {len(full_test_df) + len(full_training_df)}')

    # save training and test sets
    full_training_df.to_hdf(f"{output_dir}raw_full_training_df.h5",key="df",mode="w")
    full_test_df.to_hdf(f"{output_dir}raw_full_test_df.h5",key="df",mode="w")

if __name__ == "__main__":
    main()
