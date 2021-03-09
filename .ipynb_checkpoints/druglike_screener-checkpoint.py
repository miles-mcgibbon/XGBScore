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

def get_stats_from_MOAD_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    num_residues = len(mol.residues)
    coords = mol.coords
    if coords.size == 0:
        molecule_check = False
        num_ters = 1
        mw = 1000
        num_rotors = 14
    else:
        coords_df = pd.DataFrame({'X': coords[:, 0], 'Y': coords[:, 1], 'Z': coords[:, 2]})
        diffs = list()
        for col in coords_df:
            coords_df[f'd{col}'] = coords_df[col].diff(-1).abs()
            for diff in list(coords_df[f'd{col}']):
                diffs.append(diff)
        bond_breaks = len([d for d in diffs if d > 2.5])
        mol_text = open(filepath, 'r').read()
        num_ters = mol_text.count('TER')
        molecule_check = True if (bond_breaks/num_ters) <= 1 else False
        mw = mol.molwt/num_ters
        num_rotors = mol.num_rotors/num_ters
    return mw, num_rotors, num_residues, num_ters, molecule_check

def get_stats_from_PDB_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    num_residues = len(mol.residues)
    mw = mol.molwt/num_residues
    num_rotors = mol.num_rotors/num_residues
    num_ters = 1
    molecule_check = True
    return mw, num_rotors, num_residues, num_ters, molecule_check

def druglike_filter(mw, num_rotors, num_residues, num_ters, molecule_check):
    if mw < 510 and num_rotors <= 13:
        if num_residues == 1:
            return True
        elif molecule_check:
            return True
        else:
            return False
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

print(passes)


with tqdm(total=len(pass_folders)) as pbar:
    for filepath in pass_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = [filename for filename in os.listdir(filepath) if 'ligand.mol2' in filename or 'ligand.pdb' in filename][0]
        shutil.copytree(filepath, f'{destination_path}{foldername}')
        pbar.update(1)
