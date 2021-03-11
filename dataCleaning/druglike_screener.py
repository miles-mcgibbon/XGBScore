import os
from biopandas.mol2 import PandasMol2
from molmass import Formula
import oddt
from tqdm import tqdm
import pandas as pd
import shutil
import sys


def get_stats_from_file(filepath):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    file_format = str(filename.split('.')[1])
    mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
    mw = mol.molwt
    num_rotors = mol.num_rotors
    return mw, num_rotors

def druglike_filter(mw, num_rotors, mw_thresh, num_rotors_thresh):
    if mw <= mw_thresh and num_rotors <= num_rotors_thresh:
        return True
    else:
        return False

def parse_args(args):

    files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-des') + 1]

    mw_thresh = float(args[args.index('-mw') + 1])

    num_rotors_thresh = float(args[args.index('-rot') + 1])

    return files_path, destination_path, mw_thresh, num_rotors_thresh

def main():

    passes = 0

    fails = 0

    files_path, destination_path, mw_thresh, num_rotors_thresh = parse_args(sys.argv)

    structure_folders = [(files_path + folder) for folder in os.listdir(files_path)]

    with tqdm(total=len(structure_folders)) as pbar:
        for filepath in structure_folders:
            foldername = filepath.split('/')
            foldername = foldername[len(foldername) - 1]
            if foldername == 'index' or foldername == 'readme':
                pass
            else:
                ligand_file = [f'{filepath}/{filename}' for filename in os.listdir(filepath) if 'ligand' in filename][0]
                ligand_stats = get_stats_from_file(ligand_file)
                if druglike_filter(*ligand_stats, mw_thresh, num_rotors_thresh):
                    passes = passes + 1
                    shutil.copytree(filepath, f'{destination_path}{foldername}')
                else:
                    fails = fails + 1
            pbar.update(1)

    print(f'Filtering complete:')
    print(f'- {len(structure_folders)} total supplied ligands')
    print(f'- {passes} druglike molecules')
    print(f'- {fails} non-druglike molecules')
    print(f'- Pass percentage: {round(((passes/len(structure_folders))*100), 3)}')
    print(f'- Fail percentage: {round(((fails/len(structure_folders))*100), 3)}')

if __name__ == '__main__':
    main()