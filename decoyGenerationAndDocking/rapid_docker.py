'''
Creates pdb and pdbqt copies of DeepCoy produced .smi decoys
and then docks these pdbqt decoys into their respective
active crystal receptor.pdbqt file using GWOVina
'''

import os
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
import pprint
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from tqdm import tqdm
import shutil
import sys

def make_decoy_dict(decoy_file): # makes a dictionary with active smiles as keys and list of decoy smiles as values

    # define the dictionary for populating
    decoy_dict = dict()

    # read in the DeepCoy output file
    decoy_df = pd.read_csv(decoy_file, sep=" ", names=['Actives', 'Decoys'])

    # add the active and list of decoys to decoy_dict
    for active in decoy_df.Actives.unique():
        active_df = decoy_df.loc[decoy_df.Actives == active]
        decoy_dict[list(active_df['Actives'])[0]] = list(active_df['Decoys'])

    return decoy_dict


def make_pdbs_from_smiles(ref_df, decoy_files, decoy_pdbs): # make and save pdb copies of decoy smiles files

    # load reference .csv produced by smile_generator.py
    pdb_code_dict = dict(zip(list(ref_df['SMILE']), list(ref_df['PDB_CODE'])))

    # for each batch of decoys make pdb copies
    for decoy_file in decoy_files:
        decoy_dict = make_decoy_dict(decoy_file)
        with tqdm(total=len(decoy_dict.keys())) as pbar:
            for active in decoy_dict.keys():
                decoy_rank = 1
                pdb_code = pdb_code_dict.get(active)
                decoys = decoy_dict.get(active)
                for decoy in decoys:

                    # give the smile atoms coordinates and temporary hydrogens for addition of 3D structure
                    ligand = Chem.MolFromSmiles(decoy)
                    h_ligand = Chem.AddHs(ligand)
                    AllChem.EmbedMolecule(h_ligand,randomSeed=0xf00d)

                    # remove the temporary hydrogens from now 3D molecule
                    ligand = Chem.RemoveHs(h_ligand)

                    # save the produced pdb file
                    Chem.MolToPDBFile(ligand, f'{decoy_pdbs}{pdb_code}_decoy_{decoy_rank}.pdb')
                    decoy_rank = decoy_rank + 1
                pbar.update(1)


def autodock_convert(destination_path, prep_ligand_command, ligand_filepath): # converts files from .pdb format to .pdbqt format using AutoDockTools

    # get filename from filepath
    ligand_file = ligand_filepath.split('/')[len(ligand_filepath.split('/')) - 1]
    print('Preparing ligand...')

    # use AutoDockTools CLI to convert pdb file adding gasteiger charges and polar hydrogens
    output = subprocess.check_output(f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{ligand_file}qt -U nphs', shell=True, stderr=subprocess.STDOUT)


def make_pdbqts_from_pdbs(decoy_pdbs, decoy_pdbqts, prep_ligand_command): # make and save pdbqt copies of decoy pdb files

    # iterate over pdb decoy files and convert them
    decoy_pdb_files = [f'{decoy_pdbs}{file}' for file in os.listdir(decoy_pdbs)]
    with tqdm(total=len(decoy_pdb_files)) as pbar:
        for file in decoy_pdb_files:
            autodock_convert(decoy_pdbqts, prep_ligand_command, file)
            pbar.update(1)


def get_coordinates(file, size): # find the center x, y, z coordinates of the active crystal ligand and construct a binding site cuboid around this center point

    # read the pdb file as a dataframe
    pmol = PandasPdb().read_pdb(file)

    # find the minimum and maximum x, y, z coordinates for atoms in crystal ligand file
    x_min, x_max = min(pmol.df['ATOM'].x_coord), max(pmol.df['ATOM'].x_coord)
    y_min, y_max = min(pmol.df['ATOM'].y_coord), max(pmol.df['ATOM'].y_coord)
    z_min, z_max = min(pmol.df['ATOM'].z_coord), max(pmol.df['ATOM'].z_coord)

    # find center point for each axis
    x_center = x_min + abs(x_min-x_max)/2
    y_center = y_min + abs(y_min-y_max)/2
    z_center = z_min + abs(z_min-z_max)/2

    # calculate lengths for each axis based on difference between min and max values multiplied by user defined size variable
    x_range = (abs(x_min-x_max)+size)
    y_range = (abs(y_min-y_max)+size)
    z_range = (abs(z_min-z_max)+size)

    return x_center, y_center, z_center, x_range, y_range, z_range


def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z): # dock the decoy pdbqt to the receptor pdbqt using GWOVina CLI
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}' \
              '--exhaustiveness=32 --num_wolves=40 --num_modes=5 --energy_range=4')


def dock_all_decoys(decoy_pdbqts, pdbqt_files, docker_command, docked_decoys): # batch dock all decoy.pdbqt files in 'decoy_pdbqts' folder

    # iterate over pdbqt decoy files and dock
    decoy_pdbqt_files = [f'{decoy_pdbqts}{file}' for file in os.listdir(decoy_pdbqts)]
    with tqdm(total=len(decoy_pdbqt_files)) as pbar:
        for filepath, filename in zip(decoy_pdbqt_files, os.listdir(decoy_pdbqts)):

            # define file variables
            pdb_code = filename.split('_')[0][1:]
            foldername = filename.split('.')[0]
            receptor_file = f'{pdbqt_files}{pdb_code}/{pdb_code}_receptor.pdbqt'
            example_crystal_ligand = f'{pdbqt_files}{pdb_code}/{pdb_code}_ligand.pdbqt'

            # dock the decoy to the active crystal receptor.pdbqt
            dock_file(docker_command, receptor_file, filepath, *get_coordinates(example_crystal_ligand, 1))

            # make destination folder and transfer AutoDockTools output decoy.pdbqt to destination folder
            os.mkdir(f'{docked_decoys}{foldername}')
            shutil.copyfile(f'{decoy_pdbqts}{foldername}_out.pdbqt',f'{docked_decoys}{foldername}/{foldername}_ligand.pdbqt')
            shutil.copyfile(receptor_file, f'{docked_decoys}{foldername}/{foldername}_receptor.pdbqt')
            os.remove(f'{decoy_pdbqts}{foldername}_out.pdbqt')
            pbar.update(1)


def parse_args(args): # parse CLI user inputs

    docker_command = args[args.index('-dock') + 1]

    ref_df = pd.read_csv(args[args.index('-ref') + 1])

    pdbqt_files = args[args.index('-pdbqts') + 1]

    decoys = args[args.index('-decoys') + 1]

    docked_decoys = args[args.index('-des') + 1]

    decoy_files = [f'{decoys}{file}' for file in os.listdir(decoys)]

    prep_ligand_command = args[args.index('-prep_lig') + 1]

    return docker_command, ref_df, pdbqt_files, decoys, docked_decoys, decoy_files, prep_ligand_command


def main(): # run script using CLI

    docker_command, ref_df, pdbqt_files, decoys, docked_decoys, decoy_files, prep_ligand_command = parse_args(sys.argv)

    # make output folders if none exist
    decoy_pdbs = os.path.join(os.path.dirname(__file__), 'decoy_pdbs','')
    decoy_pdbqts = os.path.join(os.path.dirname(__file__), 'decoy_pdbqts','')
    if not os.path.isdir(decoy_pdbs):
        os.mkdir(decoy_pdbs)
    if not os.path.isdir(decoy_pdbqts):
        os.mkdir(decoy_pdbqts)

    # first convert decoys from smiles to pdbqt
    make_pdbs_from_smiles(ref_df, decoy_files, decoy_pdbs)
    make_pdbqts_from_pdbs(decoy_pdbs, decoy_pdbqts, prep_ligand_command)

    # dock the decoy pdbqts to the respetive crystal receptor.pdbqt
    dock_all_decoys(decoy_pdbqts, pdbqt_files, docker_command, docked_decoys)

if __name__ == '__main__':
    main()
