import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt

pdb_files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/PDBBind/'
moad_files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Binding_MOAD/'

pdb_structure_folders = [(pdb_files_path + folder) for folder in os.listdir(pdb_files_path)]
moad_structure_folders = [(moad_files_path + folder) for folder in os.listdir(moad_files_path)]

example_grouped_ligand = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/1ui0_1/1ui0_1_ligand.pdb'
example_split_ligand = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/6knh_1/6knh_1_ligand.pdb'

def get_stats_from_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    num_residues = len(mol.residues)
    mol_text = open(filepath, 'r').read()
    num_ters = mol_text.count('TER')
    return mol.molwt, mol.num_rotors, mol.smiles, num_residues, num_ters

def druglike_filter(mw, num_rotors, smiles, num_residues, num_ters):
    if mw < 510 and num_rotors <= 13:
        if num_residues == 1:
            return True
        elif num_residues == num_ters:
            return True
        else:
            return False
    else:
        return False

print(druglike_filter(*get_stats_from_file(example_grouped_ligand, 'pdb')))
print(druglike_filter(*get_stats_from_file(example_split_ligand, 'pdb')))
