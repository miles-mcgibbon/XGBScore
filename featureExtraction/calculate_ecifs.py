'''
This script makes an sdf copy of a pdb crystallographic
ligand then supplies this sdf ligand and its _receptor
pdb file to efic_utils function to calculate ecif
interactions between ligand and receptor atoms
'''

from ecif_utils import *
import oddt
import pandas as pd
import sys
from tqdm import tqdm

def make_sdf_copy(pdbqt_path, pdb_path, pdb_code): # uses openbabel to convert pdbqt files to sdf format

    # define filepaths
    pdbqt_sdf_ligand = f'{pdbqt_path}{pdb_code}/{pdb_code}_ligand.pdbqt'

    if not os.path.isdir(f'{pdb_path}{pdb_code}'):
        os.mkdir(f'{pdb_path}{pdb_code}')

    # convert and save in pdb copies location
    for mol in oddt.toolkits.ob.readfile('pdbqt', pdbqt_sdf_ligand):
        mol.write(format='pdb', filename=f'{pdb_path}{pdb_code}/{pdb_code}_interim_ligand.pdb', overwrite=True, opt=None, size=None)
    for mol in oddt.toolkits.ob.readfile('pdb',f'{pdb_path}{pdb_code}/{pdb_code}_interim_ligand.pdb'):
        mol.addh(only_polar=False)
        mol.write(format='sdf', filename=f'{pdb_path}{pdb_code}/{pdb_code}_sdf_ligand.sdf', overwrite=True, opt=None, size=None)



def calculate_ecif_data_entry(pdb_path, pdb_code): # calculates ecif data for specific structure and returns data as series for appending to df

    sdf_ligand = f'{pdb_path}{pdb_code}/{pdb_code}_sdf_ligand.sdf'

    clean_pdb_code = pdb_code[:4]

    pdb_receptor = f'{pdb_path}{clean_pdb_code}/{clean_pdb_code}_receptor.pdb'

    return GetECIF(pdb_receptor, sdf_ligand, distance_cutoff=6.0)

def process_structure(pdbqt_path, pdb_path, pdb_code): # multithreaded approach to making sdfs and calculating ecif data

    # make sdf copy of pdbqt ligand for ecif calculations
    make_sdf_copy(pdbqt_path, pdb_path, pdb_code)

    # calculate ecif data
    ecif_data = calculate_ecif_data_entry(pdb_path, pdb_code)

    return ecif_data

def parse_args(args): # parse CLI user inputs

    pdbqt_path = args[args.index('-pdbqts') + 1]

    pdb_path = args[args.index('-pdbs') + 1]

    save_loc = args[args.index('-out') + 1]

    return pdbqt_path, pdb_path, save_loc


def main(): # run script using CLI

    pdbqt_path, pdb_path, save_loc = parse_args(sys.argv)

    # make empty csv for populating with ecif data
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]
    columns=['PDBCode'] + ECIFHeaders
    columns = ','.join(columns)
    with open(save_loc,'a+') as master_df:
        master_df.write(f'{columns}\n')
        master_df.close()

    with tqdm(total=len(os.listdir(pdbqt_path))) as pbar:
        for pdb_code in os.listdir(pdbqt_path):
            try:
                # calculate ecifs
                ecif_data = process_structure(pdbqt_path, pdb_path, pdb_code)

                # write out row of data to df
                df_row = [pdb_code] + ecif_data
                df_row = [str(num) for num in df_row]
                df_row = ','.join(df_row)
                with open(save_loc,'a+') as master_df:
                    master_df.write(f'{df_row}\n')
                    master_df.close()
            except:
                pass
            pbar.update(1)

if __name__ == '__main__':
    main()
