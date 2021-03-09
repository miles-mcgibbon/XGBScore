import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt
from tqdm import tqdm
import pandas as pd
import shutil

files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Full/'

structure_folders = [(files_path + folder) for folder in os.listdir(files_path)]

destination_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Full_Druglike/'

def get_stats_from_file(filepath):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    file_format = str(filename.split('.')[1])
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    mw = mol.molwt
    num_rotors = mol.num_rotors
    return mw, num_rotors

def druglike_filter(mw, num_rotors):
    if mw < 510 and num_rotors <= 13:
        return True
    else:
        return False

passes = 0

fails = 0

with tqdm(total=len(structure_folders)) as pbar:
    for filepath in structure_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        if foldername == 'index' or foldername == 'readme':
            pass
        else:
            ligand_file = [f'{filepath}/{filename}' for filename in os.listdir(filepath) if 'ligand' in filename][0]
            ligand_stats = get_stats_from_file(ligand_file, 'mol2')
            if druglike_filter(*ligand_stats):
                passes = passes + 1
                shutil.copytree(filepath, f'{destination_path}{foldername}')
            else:
                fails = fails + 1
        pbar.update(1)
