'''
Fetches the pdb versions of all the structures in the pdbqt final
dataset for conversion to smiles, as pdbqts contain artificially
added polar hydrogens
'''

import os
import oddt
import subprocess
from tqdm import tqdm
import shutil
import sys

def copy_pdb_files(original_pdb_files_path, destination_path, pdbqt_path_folders): # simple function to copy pdb files to destination from source

    # prepare destination folder
    try:
        os.mkdir(destination_path)
    except FileExistsError:
        pass

    # copy folder
    with tqdm(total=len(pdbqt_path_folders)) as pbar:
        for folder in pdbqt_path_folders:
            folder_name = folder.split('/')
            folder_name = folder_name[len(folder_name) - 1]
            shutil.copytree(f'{original_pdb_files_path}{folder_name}', f'{destination_path}{folder_name}')
            pbar.update(1)

def parse_args(args): # parse CLI user inputs

    original_pdb_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-des') + 1]

    pdbqt_path = args[args.index('-ref') + 1]

    return original_pdb_files_path, destination_path, pdbqt_path

def main(): # run script using CLI

    original_pdb_files_path, destination_path, pdbqt_path = parse_args(sys.argv)

    pdbqt_path_folders = [(pdbqt_path + folder) for folder in os.listdir(pdbqt_path)]

    copy_pdb_files(original_pdb_files_path, destination_path, pdbqt_path_folders)

if __name__ == '__main__':
    main()
