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

    # convert and save in pdb copies location
    for mol in oddt.toolkits.ob.readfile('pdbqt', pdbqt_sdf_ligand):
        mol.write(format='sdf', filename=f'{pdb_path}{pdb_code}/{pdb_code}_sdf_ligand.sdf', overwrite=True, opt=None, size=None)


def calculate_ecif_data_entry(pdb_path, pdb_code): # calculates ecif data for specific structure and returns data as series for appending to df
    sdf_ligand = f'{pdb_path}{pdb_code}/{pdb_code}_sdf_ligand.sdf'
    pdb_receptor = f'{pdb_path}{pdb_code}/{pdb_code}_receptor.pdb'

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

    # make empty dataframe for populating with ecif data
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]
    master_df = pd.DataFrame(columns=(['PDBCode'] + ECIFHeaders))


    with tqdm(total=len(os.listdir(pdbqt_path))) as pbar:
        for pdb_code in os.listdir(pdbqt_path):

            # calculate ecifs
            ecif_data = process_structure(pdbqt_path, pdb_path, pdb_code)

            # add row of data to df
            df_row = [pdb_code] + ecif_data
            master_df.loc[len(master_df)] = df_row
            pbar.update(1)

    master_df.to_csv(save_loc, index=False)

if __name__ == '__main__':
    main()
