import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt

example_crystal_ligand = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/crystal_ligand.mol2'

def get_mw_from_file(filepath, file_format):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    mol2 = next(oddt.toolkits.ob.readfile(file_format, filepath))
    return mol2.molwt

print(get_mw_from_file(example_crystal_ligand, 'mol2'))
