'''
This script gets the binding data values from the database binding
data file that the structure was sourced from, produces a csv of
structures with binding data that can be classified with a user
defined Kd threshold, and offers to make a copy of the dataset
only containing structures with classified binding data
'''

import pandas as pd
import sys
import numpy as np
from io import StringIO
import os
import shutil
from tqdm import tqdm
from warnings import filterwarnings
filterwarnings("ignore")

pd.set_option('display.max_columns', None)

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

def unit_convert(df): # standardise all binding data to uM units

    print('Converting units...\n')
    df['Binding_Data_in_uM'] = df['BindingValue'].apply(lambda x:  x[0:] if x.find("(") == -1 else x[0:x.find("(")]).astype(float)

    # slice dataframes into units and convert to uM
    kadf = df[df['BindingDataType'] == 'Ka']
    kadf['Binding_Data_in_uM'] = (1/kadf['Binding_Data_in_uM']) * 1e6

    micro = df[df['BindingUnits'] == 'uM']

    nano = df[df['BindingUnits'] == 'nM']
    nano['Binding_Data_in_uM'] = nano['Binding_Data_in_uM'] * 1e-3

    pico = df[df['BindingUnits'] == 'pM']
    pico['Binding_Data_in_uM'] = pico['Binding_Data_in_uM']* 1e-6

    femto = df[df['BindingUnits'] == 'fM']
    femto['Binding_Data_in_uM'] = femto['Binding_Data_in_uM']* 1e-9

    molar = df[df['BindingUnits'] == 'M']
    molar['Binding_Data_in_uM'] = molar['Binding_Data_in_uM']* 1e6

    milli = df[df['BindingUnits'] == 'mM']
    milli['Binding_Data_in_uM'] = milli['Binding_Data_in_uM'] * 1e3

    # concat converted uM dataframes to master dataframe
    all = pd.concat([kadf,micro,nano,pico,femto,molar,milli]).sort_index(ascending=True)

    return all

def cc50_check(df): # remove values stored as CC50 as these cannot be logically converted to Kd

    print('\nChecking for CC50 values...')

    # isolate CC50 values in the df
    cc = df[df['BindingDataType'].str.contains('CC50' or 'cc50' or 'Cc50')]

    # set up list for populating
    ccPDB = list()

    # populate CC50 list with pdb codes of CC50 structures
    if len(cc) > 0:
        print('File contains CC50 values')
        print('Removing...')
        df = df[df['BindingDataType'] != 'CC50']
        df = df[df['BindingDataType'] != 'cc50']
        df = df[df['BindingDataType'] != 'Cc50'].reset_index()
        print('All CC50 values:\n')
        print(cc)
        ccPDB = list(cc['PDBCode'])

    return df, ccPDB

def ic50_check(df, threshold): # remove IC50 values outside of user defined threshold

    print(f'Filtering out IC50 values between {threshold} uM and {threshold*10} uM...')

    # keep structures below or above the threshold
    ic = df[df['BindingDataType'].str.contains('50')]
    ic_below = ic[ic['Binding_Data_in_uM'].astype(float) <= threshold]
    ic_above = ic[ic['Binding_Data_in_uM'].astype(float) >= threshold*10]
    ic_kept = pd.concat([ic_above,ic_below])

    # define pdb codes of structures to keep
    kept = list(ic_kept['PDBCode'])

    # define pdb codes of structures to remove
    icPDB = list(ic['PDBCode'])
    removed = [structure for structure in icPDB if structure not in kept]

    # add kept IC50 structures back to master dataframe
    df = df[~df['BindingDataType'].str.contains('50')]
    df = pd.concat([df,ic_kept]).sort_index(ascending=True)

    return df, removed

def parse_args(args): # parse CLI user inputs

    iridium_data = args[args.index('-irid') + 1]

    pdb_bind_data = args[args.index('-pdb') + 1]

    binding_moad_data = args[args.index('-moad') + 1]

    pdbqt_files_path = args[args.index('-pdbqt') + 1]

    threshold = float(args[args.index('-threshold') + 1])

    output_dir = args[args.index('-out') +1]

    return iridium_data, pdb_bind_data, binding_moad_data, pdbqt_files_path, threshold, output_dir


def main(): # run script using CLI

    iridium_data, pdb_bind_data, binding_moad_data, pdbqt_files_path, threshold, output_dir = parse_args(sys.argv)

    # collate raw binding data
    master_df = load_and_concat_all_data(iridium_data, pdb_bind_data, binding_moad_data)

    # build dataframe of structures with any binding data present
    df, missing_structures = search_for_binding_data(master_df, pdbqt_files_path)

    # save summary csv for reference to script location and print errors
    df.to_csv(f'{output_dir}bindingDataSummary.csv', index=False)

    # remove CC50 structures
    df['BindingValue'] = df['BindingValue'].astype(str)
    df, ccPDB = cc50_check(df)

    # convert all binding data to uM
    df = unit_convert(df)

    # remove IC50 structures outside of the user defined threshold values
    df, removed = ic50_check(df, threshold)

    # add classification label to binding data based on user threshold
    df['Label'] = np.where(df['Binding_Data_in_uM'] < threshold, 1, 0)

    # save useful binding data with kd label
    df.pop('index')
    df.to_csv(f'{output_dir}convertedBindingData.csv', index=False)

    # update user and give summary of data cleaning steps
    if len(removed) == 0 and len(ccPDB) == 0:
        print('No CC50 values found')
        print('No IC50 values removed.')
    else:
        print(f'{len(ccPDB)} structure(s) with CC50 values:')
        print(ccPDB)
        print(f'{len(removed)} structure(s) do not meet IC50 value requirements:')
        print(removed)
        print(f'{len(missing_structures)} structure(s) have no binding data present:')
        print(missing_structures)

        # option to make copy of dataset with actives that have classified binding data only
        answer = input('Make copy of dataset without these structures? (y/n)')
        if answer == 'y':
            os.mkdir(f'{pdbqt_files_path[:len(pdbqt_files_path) - 1]}_filtered_binding_data')
            print('Copying files...')
            with tqdm(total=len(df)) as pbar:
                for structure in list(df['PDBCode']):
                    shutil.copytree(f'{pdbqt_files_path}{structure}', f'{pdbqt_files_path[:len(pdbqt_files_path) - 1]}_filtered_binding_data/{structure}')
                    pbar.update(1)

if __name__ == '__main__':
    main()
