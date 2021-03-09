import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt
from tqdm import tqdm
import pandas as pd
import shutil

pdb_files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/PDBBind/'
moad_files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Binding_MOAD/'

pdb_structure_folders = [(pdb_files_path + folder) for folder in os.listdir(pdb_files_path)]
moad_structure_folders = [(moad_files_path + folder) for folder in os.listdir(moad_files_path)]


destination_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Full_Druglike/'

def get_stats_from_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
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

pass_folders = list()

with tqdm(total=len(moad_structure_folders)) as pbar:
    for filepath in moad_structure_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = filepath + '/' + foldername + '_ligand.pdb'
        ligand_stats = get_stats_from_MOAD_file(ligand_file, 'pdb')
        if druglike_filter(*ligand_stats):
            passes = passes + 1
            pass_folders.append(filepath)
        else:
            fails = fails + 1
        pbar.update(1)

with tqdm(total=len(pdb_structure_folders)) as pbar:
    for filepath in pdb_structure_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        if foldername == 'index' or foldername == 'readme':
            fails = fails + 1
        else:
            ligand_file = filepath + '/' + foldername + '_ligand.mol2'
            ligand_stats = get_stats_from_PDB_file(ligand_file, 'mol2')
            if druglike_filter(*ligand_stats):
                passes = passes + 1
                pass_folders.append(filepath)
            else:
                fails = fails + 1
        pbar.update(1)

with tqdm(total=len(pass_folders)) as pbar:
    for filepath in pass_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = [filename for filename in os.listdir(filepath) if 'ligand.mol2' in filename or 'ligand.pdb' in filename][0]
        shutil.copytree(filepath, f'{destination_path}{foldername}')
        pbar.update(1)
