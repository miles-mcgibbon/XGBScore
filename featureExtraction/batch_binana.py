'''
This script batch processes receptor and ligand pdbqt files and
passes them to binana.py to produce text files containing
descriptors of the protein-ligand interaction
'''

import os
import shutil
import subprocess
import sys
from tqdm import tqdm
import concurrent.futures
from subprocess import PIPE

fatal_error_list = list()

binana_path = None
destination_path = None

def run_binana(folder_name, ligand_path, receptor_path, ligand): # use binana.py with command line arguments as the ligand and receptor pdbqts and save text file to destination path

    # run binana on ligand receptor pdbqt pair
    try:
        output = subprocess.check_output(f'python {binana_path} -receptor {receptor_path} -ligand {ligand_path} > {destination_path}{ligand}.txt', shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:

        # record error and save to fatal_error_list
        fatal_error_list.append(folder_name)
        if os.path.isfile('fatal_error_list.txt'):
            os.remove('fatal_error_list.txt')
        with open('fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list) + '\n' + str(e))
            error_list.close()
        print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

def extract_binana_features(folder_path): # wrapper function for passing structure pdbqt files to binana.py

    # define filepaths
    folder_name = folder_path.split('/')
    folder_name = folder_name[len(folder_name) - 1]

    # define ligand filepath and ligand variables
    ligand_file = [filename for filename in os.listdir(folder_path) if 'ligand.pdbqt' in filename][0]
    ligand_path = folder_path + '/' + ligand_file
    ligand = ligand_file.replace('.pdbqt','')

    # define receptor filepath and receptor variables
    receptor_file = [filename for filename in os.listdir(folder_path) if 'receptor.pdbqt' in filename][0]
    receptor_path = folder_path + '/' + receptor_file

    # calculate binana features
    run_binana(folder_name, ligand_path, receptor_path, ligand)

def batch_extract_binana_features(pdbqt_files,pdbqt_files_path ,threads): # multithread the processing of structures by binana.py
    threads = min(threads, len(pdbqt_files))
    print(f'Running BINANA on {pdbqt_files_path}...')
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(extract_binana_features, pdbqt_files), total=len(pdbqt_files)))


def move_problematic_structures(fatal_error_list, pdbqt_files_path, problem_path): # copy structures that threw exceptions to folder for debugging

    # remove duplicates from fatal error list
    fatal_error_list = list(set(list(fatal_error_list)))
    print(fatal_error_list)

    # copy folders to user defined problem directory
    with tqdm(total=len(fatal_error_list)) as pbar:
        for folder in fatal_error_list:
            shutil.move(f'{pdbqt_files_path}{folder}', f'{problem_path}{folder}')
            #shutil.copy(f'{destination_path}{folder}.txt',f'{problem_path}{folder}/{folder}.txt')
            #os.remove(f'{destination_path}{folder}.txt')
            pbar.update(1)


def parse_args(args): # parse CLI user inputs

    binana_path = args[args.index('-binana') + 1]

    pdbqt_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-suc') + 1]

    problem_path = args[args.index('-prob') + 1]

    threads = int(args[args.index('-threads') + 1])

    return binana_path, pdbqt_files_path, destination_path, problem_path, threads

def main(): # run script using CLI

    global destination_path
    global binana_path

    binana_path, pdbqt_files_path, destination_path, problem_path, threads = parse_args(sys.argv)

    # make absolute paths for structure folders
    pdbqt_folders = [(pdbqt_files_path + folder) for folder in os.listdir(pdbqt_files_path)]

    # calculate binana features
    batch_extract_binana_features(pdbqt_folders, pdbqt_files_path, threads)

    move_problematic_structures(fatal_error_list, pdbqt_files_path, problem_path)

if __name__ == '__main__':
    main()
