import os
from tqdm import tqdm
import pypdb
import pprint
import pandas as pd
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def convert_file_to_smiles(filepath, filetype, phosphorous_only):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    if filetype == 'pdb':
        mol = Chem.MolFromPDBFile(filepath)
        smi = Chem.MolToSmiles(mol)
    elif filetype == 'mol2':
        mol = Chem.MolFromMol2File(filepath)
        smi = Chem.MolToSmiles(mol)
    if len(smi.strip().split()) == 0:
        raise BadLigandException
    if phosphorous_only:
        if 'P' in str(smi):
            return smi
        elif 'p' in str(smi):
            return smi
        else:
            return None
    else:
        if 'P' in str(smi):
            return None
        elif 'p' in str(smi):
            return None
        else:
            return smi

def load_active_crystal_ligands(pdb_folders, backup_folder, phosphorous_only):
    crystal_ligands_dict = dict()
    errors = list()
    for filepath in pdb_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = [(filepath + '/' + filename) for filename in os.listdir(filepath) if 'ligand.pdb' in filename][0]
        try:
            ligand_smile = convert_file_to_smiles(ligand_file, 'pdb', phosphorous_only).strip().lstrip()
            if ligand_smile is not None:
                crystal_ligands_dict[foldername] = ligand_smile
        except:
            try:
                ligand_smile = convert_file_to_smiles(f'{backup_folder}{foldername}/{foldername}_ligand.mol2', 'mol2', phosphorous_only).strip().lstrip()
                if ligand_smile is not None:
                    crystal_ligands_dict[foldername] = ligand_smile
            except Exception as e:
                print(str(e))
                errors.append(foldername)
    return crystal_ligands_dict, errors


def parse_args(args):

    pdb_druglike_files_path = args[args.index('-loc') + 1]

    backup_path = args[args.index('-backup') + 1]

    split_files = args[args.index('-split') + 1]

    if '-phos' in args:
        phosphorous_only = True
    else:
        phosphorous_only = False

    return pdb_druglike_files_path, backup_path, int(split_files), phosphorous_only

def main():

    pdb_druglike_files_path, backup_path, split_files, phosphorous_only = parse_args(sys.argv)

    pdb_folders = [(pdb_druglike_files_path + folder) for folder in os.listdir(pdb_druglike_files_path)]

    actives_dict, errors = load_active_crystal_ligands(pdb_folders, backup_path, phosphorous_only)

    print(f'Completed with {len(errors)} skipped files:')
    print(errors)

    # make file for DeepCoy
    if phosphorous_only:
        # make reference .csv file
        df = pd.DataFrame({'PDB_CODE':actives_dict.keys(), 'SMILE':actives_dict.values()})
        df.to_csv('phosphorousSmileReferenceSheet.csv')

        smiles = list(actives_dict.values())
        print(len(smiles))
        split_files = int(len(smiles)/split_files)
        split_smiles = [smiles[i:i + split_files] for i in range(0, len(smiles), split_files)]
        smile_markers = [f'{i}-{i+split_files}' for i in range(0, len(smiles), split_files)]
        for smile_list, marker in zip(split_smiles, smile_markers):
            outfile = open(f"p_actives_{marker}.smi", "w+")
            outfile.write("\n".join(str(i) for i in smile_list))
            outfile.close()

    else:
        smiles = list(actives_dict.values())
        print(len(smiles))
        split_files = int(len(smiles)/split_files)
        split_smiles = [smiles[i:i + split_files] for i in range(0, len(smiles), split_files)]
        smile_markers = [f'{i}-{i+split_files}' for i in range(0, len(smiles), split_files)]
        for smile_list, marker in zip(split_smiles, smile_markers):
            outfile = open(f"np_actives_{marker}.smi", "w+")
            outfile.write("\n".join(str(i) for i in smile_list))
            outfile.close()

        # make reference .csv file
        df = pd.DataFrame({'PDB_CODE':actives_dict.keys(), 'SMILE':actives_dict.values()})
        df.to_csv('noPhosphorousSmileReferenceSheet.csv')


if __name__ == '__main__':
    main()
