from oddt.scoring.descriptors.binana import binana_descriptor
import oddt
import os
from tqdm import tqdm
import pypdb
import pprint
import pandas as pd

pdbqt_druglike_files_path = '/home/milesm/Dissertation/Data/PDBQT/Non_redundant/Druglike/'

pdbqt_folders = [(pdbqt_druglike_files_path + folder) for folder in os.listdir(pdbqt_druglike_files_path)]

def convert_file_to_smiles(filepath, filetype):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    molecule = next(oddt.toolkits.ob.readfile(filetype, filepath))
    return molecule.smiles

def try_binana(pdbqt_folders):
    for filepath in pdbqt_folders:
        if '2wxd' in filepath:
            foldername = filepath.split('/')
            foldername = foldername[len(foldername) - 1]
            print(foldername)
            ligand_file = [(filepath + '/' + filename) for filename in os.listdir(filepath) if 'ligand.pdbqt' in filename][0]
            receptor_file = [(filepath + '/' + filename) for filename in os.listdir(filepath) if 'receptor.pdbqt' in filename][0]
            ligand = oddt.toolkits.ob.readfile('pdbqt', ligand_file)
            receptor = next(oddt.toolkits.ob.readfile('pdbqt', receptor_file))
            print(binana_descriptor(receptor).build(ligand))
            input('')

try_binana(pdbqt_folders)
