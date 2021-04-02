'''
This script scrapes all binding data for the pdbqt dataset
from the respective binding data files supplied with the
datasets
'''

import pandas as pd
from io import StringIO
import os
import sys
import shutil
from tqdm import tqdm

def load_moad(binding_moad_data): # load binding data from Binding MOAD dataset into dataframe

    # load binding MOAD reference dataframe with forced headers and clean up formatting
    force_columns = ['Decimal',
                 'EntryString',
                 'PDBCode',
                 'LigandID',
                 'Validity',
                 'BindingDataType',
                 'Equals',
                 'BindingValue',
                 'BindingUnits',
                 'FormulaString',
                 'FormulaStringSpill']
    binding_df = pd.read_csv(binding_moad_data, names=force_columns)
    binding_df['PDBCode'] = binding_df['PDBCode'].ffill()

    # drop rows not containing binding data
    binding_df.dropna(how='any', subset=['BindingDataType'], inplace=True)
    binding_df.drop_duplicates(subset=['PDBCode'], inplace=True)
    binding_df['Database'] = 'MOAD'

    # remove redundant columns and set formatting
    binding_df = binding_df[['PDBCode','BindingDataType','BindingValue','BindingUnits','Database']]
    binding_df['PDBCode'] = binding_df['PDBCode'].str.lower()

    return binding_df

def load_pdbbind(pdb_bind_data): # load binding data from PDBind dataset into dataframe

    # skip the first six lines with explanation of data
    with open(pdb_bind_data,'r+') as text_data:
        lines = text_data.read()
    lines = lines.split('\n')
    lines = lines[6:]

    # remove all spaces in lines
    lines = [line.replace('    ',';') for line in lines]
    lines = [line.replace('   ',';') for line in lines]
    lines = [line.replace('  ',';') for line in lines]
    lines = [line.replace(' ',';') for line in lines]
    lines = '\n'.join(lines)

    # create dataframe from text file of lines
    lines = StringIO(lines)
    df = pd.read_csv(lines, sep=';', names=['PDBCode','Resolution','Release Year','-logK','K','none','none_2','code_2','lig_name'])

    # split binding data into separate columns and add to dataframe
    binding_data = df['K']
    data_type = [data.split('=')[0] for data in binding_data]
    data_value = [data.split('=')[1][:len(data.split('=')[1]) - 2] for data in binding_data]
    data_units = [data[len(data) - 2:] for data in binding_data]
    df['BindingDataType'] = data_type
    df['BindingValue'] = data_value
    df['BindingUnits'] = data_units
    df['Database'] = 'PDBBind'

    # remove redundant columns and set formatting
    df = df[['PDBCode','BindingDataType','BindingValue','BindingUnits','Database']]
    df['PDBCode'] = df['PDBCode'].str.lower()

    return df

def load_iridium(iridium_data): # load binding data from Iridium dataset into dataframe

    # load the data
    df = pd.read_excel(iridium_data)

    # drop empty rows/merged cells and add necessary columns
    df.dropna(subset=['Binding affinity (uM)'], inplace=True)
    df.rename(columns={'Binding affinity (uM)':'BindingValue','Type of binding':'BindingDataType','Protein-Ligand':'PDBCode'}, inplace=True)
    df['Database'] = 'Iridium'
    df['BindingUnits'] = 'uM'

    # remove redundant columns and set formatting
    df = df[['PDBCode','BindingDataType','BindingValue','BindingUnits','Database']]
    df['PDBCode'] = df['PDBCode'].str.lower()

    return df

def load_and_concat_all_data(iridium_data, pdb_bind_data, binding_moad_data): # make a master dataframe of all the binding data from all three databases

    # load all the data from all databases
    moad_df = load_moad(binding_moad_data)
    pdb_df = load_pdbbind(pdb_bind_data)
    iridium_df = load_iridium(iridium_data)

    # concat these into a master dataframe
    master_df = pd.concat([moad_df, pdb_df, iridium_df])

    return master_df

def search_for_binding_data(master_df, pdbqt_files_path): # assign binding data to each structure in our dataset and return a dataframe

    # get a list of structures in dataset
    structures = os.listdir(pdbqt_files_path)

    # setup list for populating with dataframes
    dfs = list()

    # search in preference order iridium, pdbbind, moad and add structures to dfs list
    for database in ['Iridium','PDBBind','MOAD']:
        database_df = master_df.loc[master_df.Database == database]
        found_df = database_df.loc[database_df.PDBCode.isin(structures)]
        found_structures = list(found_df['PDBCode'])
        structures = [structure for structure in structures if structure not in found_structures]
        dfs.append(found_df)

    # combine the dataframes
    binding_data_final = pd.concat(dfs)

    # collect structures with no binding data
    found = list(binding_data_final['PDBCode'])
    structures = os.listdir(pdbqt_files_path)
    missing_structures = [structure for structure in structures if structure not in found]

    return binding_data_final, missing_structures

def parse_args(args): # parse CLI user inputs

    iridium_data = args[args.index('-irid') + 1]

    pdb_bind_data = args[args.index('-pdb') + 1]

    binding_moad_data = args[args.index('-moad') + 1]

    pdbqt_files_path = args[args.index('-pdbqt') + 1]

    return iridium_data, pdb_bind_data, binding_moad_data, pdbqt_files_path

def main(): # run script using CLI

    # get user args
    iridium_data, pdb_bind_data, binding_moad_data, pdbqt_files_path = parse_args(sys.argv)

    # collate raw binding data
    master_df = load_and_concat_all_data(iridium_data, pdb_bind_data, binding_moad_data)

    # match to structures in the dataset
    final_df, missing_structures = search_for_binding_data(master_df, pdbqt_files_path)

    # save output to script location and print errors
    final_df.to_csv(os.path.join(os.path.dirname(__file__), 'bindingData.csv'), index=False)
    if len(missing_structures) == 0:
        print('Binding data found for all structures')
    else:
        print(f'Unable to find binding data for {len(missing_structures)} structures:')
        print(missing_structures)
        answer = input('Make copy of dataset without these structures? (y/n)')
        if answer == 'y':
            os.mkdir(f'{pdbqt_files_path[:len(pdbqt_files_path) - 1]}_with_binding_data')
            print('Copying files...')
            with tqdm(total=len(final_df)) as pbar:
                for structure in list(final_df['PDBCode']):
                    shutil.copytree(f'{pdbqt_files_path}{structure}', f'{pdbqt_files_path[:len(pdbqt_files_path) - 1]}_with_binding_data/{structure}')
                    pbar.update(1)

if __name__ == '__main__':
    main()
