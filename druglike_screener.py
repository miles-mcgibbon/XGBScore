import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt

def get_stats_from_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    return mol.molwt, mol.num_rotors, mol.smiles

def druglike_filter(mw, num_rotors, smiles):
    if mw < 510 and num_rotors <= 13 and len(smiles.split('.') == 1):
        return True
    else:
        return False
