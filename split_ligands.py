import os
from biopandas.mol2 import *
from tqdm import tqdm
from warnings import filterwarnings
from shutil import copyfile
import requests


def fetch_missing_decoys(foldername, ligands_path, decoys_filename):

    # define the url where the decoys mol2.gz file is stored
    url = f'http://dude.docking.org/targets/{foldername}/{decoys_filename}'

    # define where to save the file
    target_path = f'{ligands_path}{foldername}/{decoys_filename}'

    # ping the url and download the file
    response = requests.get(url, verify=False)
    if response.status_code == 200:
        with open(target_path, 'wb') as file:
            file.write(response.content)
            file.close()
            print('Fetched decoys!')
    else:
        print('Error fetching decoys: Passing...')

def fetch_missing_actives(foldername, ligands_path, actives_filename):

    # define the url where the actives mol2.gz file is stored
    url = f'http://dude.docking.org/targets/{foldername}/{actives_filename}'

    # define where to save the file
    target_path = f'{ligands_path}{foldername}/{actives_filename}'

    # ping the url and download the file
    response = requests.get(url, verify=False)
    if response.status_code == 200:
        with open(target_path, 'wb') as file:
            file.write(response.content)
            file.close()
            print('Fetched actives!')
    else:
        print('Error fetching actives: Passing...')

def split_multimol_file(foldername, ligands_path, actives_filename, decoys_filename, destination_path, make_dir):


    # initialise biopandas instance
    pdmol = PandasMol2()

    # make target directories if they don't exist
    if make_dir:
        try:
            os.mkdir(f'{destination_path}{foldername}')
            os.mkdir(f'{destination_path}{foldername}/decoys')
            os.mkdir(f'{destination_path}{foldername}/actives')
        except FileExistsError:
            pass


    # split the active ligands into separate mol2 files for docking
    try:
        # try to find the actives in the existing data
        for mol2 in split_multimol2(f'{ligands_path}{foldername}/{actives_filename}'):

            pdmol.read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])

            # split the mol2 multifile and save the individual ligands
            with open(f'{destination_path}{foldername}/actives/{pdmol.code}.mol2', 'w') as file:
                file.write(pdmol.mol2_text)
                file.close()

    except FileNotFoundError:

        # fetch the ligand .gz file from dud-e and save it in the target directory
        print(f'No actives found for {foldername}, fetching from DUD-E...')
        fetch_missing_actives(foldername, ligands_path, actives_filename)

        for mol2 in split_multimol2(f'{ligands_path}{foldername}/{actives_filename}'):

            pdmol.read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])

            # split the mol2 multifile and save the individual ligands
            with open(f'{destination_path}{foldername}/actives/{pdmol.code}.mol2', 'w') as file:
                file.write(pdmol.mol2_text)
                file.close()

    # split the decoy ligands into separate mol2 files for docking
    try:
        # try to find the decoys from the existing data
        for mol2 in split_multimol2(f'{ligands_path}{foldername}/{decoys_filename}'):

            pdmol.read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])

            # split the mol2 multifile and save the individual decoys
            with open(f'{destination_path}{foldername}/decoys/{pdmol.code}.mol2', 'w') as file:
                file.write(pdmol.mol2_text)
                file.close()

    except FileNotFoundError:

        # fetch the decoy .gz file from dud-e and save it in the target directory
        print(f'No decoys found for {foldername}, fetching from DUD-E...')
        fetch_missing_decoys(foldername, ligands_path, decoys_filename)

        for mol2 in split_multimol2(f'{ligands_path}{foldername}/{decoys_filename}'):

            pdmol.read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])

            with open(f'{destination_path}{foldername}/decoys/{pdmol.code}.mol2', 'w') as file:
                file.write(pdmol.mol2_text)
                file.close()

    # copy over the receptor and crystal_ligand files to the new directory
    for file in ['crystal_ligand.mol2','receptor.pdb']:
        try:
            copyfile(f'{ligands_path}{foldername}/{file}', f'{destination_path}{foldername}/{file}')
        except FileNotFoundError:
            pass


ligands_path = '/home/milesm/Dissertation/Data/Raw/DUD-E/Extracted/'
actives_filename = 'actives_final.mol2.gz'
decoys_filename = 'decoys_final.mol2.gz'
destination_path = '/home/milesm/Dissertation/Data/Parsed/DUD-E/'

subfolders = os.listdir(ligands_path)

with tqdm(total=len(subfolders)) as pbar:
    for folder in subfolders:
        try:
            split_multimol_file(folder, ligands_path, actives_filename, decoys_filename, destination_path, True)
        except:
            print(f'Unresolvable error with {folder}')
        pbar.update(1)
