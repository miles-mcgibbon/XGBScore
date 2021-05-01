'''
Converts all active ligands in the dataset into smiles format
for decoy generation with DeepCoy, and produces a reference
csv file relating pdb codes to active smile strings
'''

import os
from tqdm import tqdm
import pandas as pd
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def convert_file_to_smiles(filepath, filetype, phosphorous_only): # use rdkit to convert ligand file to smiles format

    # get variables from filename
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]

    # convert to smiles file
    if filetype == 'pdb':
        mol = Chem.MolFromPDBFile(filepath)
        smi = Chem.MolToSmiles(mol)
    elif filetype == 'mol2':
        mol = Chem.MolFromMol2File(filepath)
        smi = Chem.MolToSmiles(mol)

    # check for empty ligands and raise exception
    if len(smi.strip().split()) == 0:
        raise BadLigandException

    if phosphorous_only:

        # return only phosphorous containing ligands
        if 'P' in str(smi):
            return smi
        elif 'p' in str(smi):
            return smi
        else:
            return None
    else:

        # return only ligands that do not contain phosphorous
        if 'P' in str(smi):
            return None
        elif 'p' in str(smi):
            return None
        else:
            return smi

def load_active_crystal_ligands(pdb_folders, backup_folder, phosphorous_only): # construct dictionary of pdb_code keys and active ligand smiles values

    # define variables for populating
    crystal_ligands_dict = dict()
    errors = list()

    # loop through pdb structure folders
    for filepath in pdb_folders:

        # get variables from filename
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = [(filepath + '/' + filename) for filename in os.listdir(filepath) if 'ligand.pdb' in filename][0]

        # try to convert pdb version to smile and add to dictionary
        try:
            ligand_smile = convert_file_to_smiles(ligand_file, 'pdb', phosphorous_only).strip().lstrip()
            if ligand_smile is not None:
                crystal_ligands_dict[foldername] = ligand_smile

        # try original mol2 file for PDBBind ligands which raise errors
        except:
            try:
                ligand_smile = convert_file_to_smiles(f'{backup_folder}{foldername}/{foldername}_ligand.mol2', 'mol2', phosphorous_only).strip().lstrip()
                if ligand_smile is not None:
                    crystal_ligands_dict[foldername] = ligand_smile

            # add ligands with unresolvable errors to error list
            except Exception as e:
                print(str(e))
                errors.append(foldername)

    return crystal_ligands_dict, errors

def parse_args(args): # parse CLI user inputs

    pdb_druglike_files_path = args[args.index('-loc') + 1]

    backup_path = args[args.index('-backup') + 1]

    split_files = args[args.index('-split') + 1]

    if '-phos' in args:
        phosphorous_only = True
    else:
        phosphorous_only = False

    return pdb_druglike_files_path, backup_path, int(split_files), phosphorous_only

def main(): # run script using CLI

    pdb_druglike_files_path, backup_path, split_files, phosphorous_only = parse_args(sys.argv)

    pdb_folders = [(pdb_druglike_files_path + folder) for folder in os.listdir(pdb_druglike_files_path)]

    # create dictionary of pdb_codes and active smiles
    actives_dict, errors = load_active_crystal_ligands(pdb_folders, backup_path, phosphorous_only)

    print(f'Completed with {len(errors)} skipped files:')
    print(errors)

    # make file for DeepCoy
    if phosphorous_only:
        # make reference .csv file
        df = pd.DataFrame({'PDB_CODE':actives_dict.keys(), 'SMILE':actives_dict.values()})
        df.to_csv('phosphorousSmileReferenceSheet.csv')

        # write out the actives smiles for DeepCoy
        smiles = list(actives_dict.values())
        print(len(smiles))
        split_files = int(len(smiles)/split_files)
        split_smiles = [smiles[i:i + split_files] for i in range(0, len(smiles), split_files)]
        smile_markers = [f'{i}-{i+split_files}' for i in range(0, len(smiles), split_files)]
        for smile_list, marker in zip(split_smiles, smile_markers):
            outfile = open(f"p_smis/p_actives_{marker}.smi", "w+")
            outfile.write("\n".join(str(i) for i in smile_list))
            outfile.close()

    else:
        # make reference .csv file
        df = pd.DataFrame({'PDB_CODE':actives_dict.keys(), 'SMILE':actives_dict.values()})
        df.to_csv('noPhosphorousSmileReferenceSheet.csv')

        # write out the actives smiles for DeepCoy
        smiles = list(actives_dict.values())
        print(len(smiles))
        split_files = int(len(smiles)/split_files)
        print(split_files)
        split_smiles = [smiles[i:i + split_files] for i in range(0, len(smiles), split_files)]
        smile_markers = [f'{i}-{i+split_files}' for i in range(0, len(smiles), split_files)]
        for smile_list, marker in zip(split_smiles, smile_markers):
            outfile = open(f"np_smis/np_actives_{marker}.smi", "w+")
            outfile.write("\n".join(str(i) for i in smile_list))
            outfile.close()

if __name__ == '__main__':
    main()
