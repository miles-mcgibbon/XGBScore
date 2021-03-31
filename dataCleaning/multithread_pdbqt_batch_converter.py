'''
Batch converts pdb receptor and ligand files to pdbqt format
across multiple cores using AutoDockTools scripts, adding
gasteiger charges and polar hydrogens
'''

import os
import oddt
import subprocess
import concurrent.futures
from tqdm import tqdm
import shutil
import sys

# define list for unresolvable errors
fatal_error_list = list()

# temporary fix for global variable error
destination_path = None
prep_ligand_command = None
prep_protein_command = None


def autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file, filetype): # converts pdb files to pdbqt

    # prepare destination folder
    try:
        os.mkdir(f'{destination_path}{folder_name}')
        print(f'Made directory {destination_path}{folder_name}')
    except FileExistsError:
        pass
    print(folder_name)

    # get filetype and alter ligand outfile
    if filetype == 'mol2':
        ligand_outfile = ligand_file.replace('mol2','pdb')
    else:
        ligand_outfile = ligand_file

    # convert ligand to pdbqt using AutoDockTools adding polar hydrogens and gasteiger charges
    print('Preparing ligand...')
    try:
        output = subprocess.check_output(f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_outfile}qt -U nphs', shell=True, stderr=subprocess.STDOUT)
    # add problematic ligands to fatal_error_list
    except:
        fatal_error_list.append(folder_name)
        if os.path.isfile('fatal_error_list.txt'):
            os.remove('fatal_error_list.txt')
        with open('fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list))
            error_list.close()
        print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

    # convert receptor to pdbqt using AutoDockTools adding polar hydrogens and gasteiger charges
    print('Preparing protein...')
    try:
        output = subprocess.check_output(f'{prep_protein_command} -r {receptor_filepath} -A hydrogens -o {destination_path}{folder_name}/{receptor_file}qt -U nphs', shell=True)

    # add problematic ligands to fatal_error_list
    except:
        fatal_error_list.append(folder_name)
        if os.path.isfile('fatal_error_list.txt'):
            os.remove('fatal_error_list.txt')
        with open('fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list))
            error_list.close()
        print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

def convert_to_pdbqt(folder): # wrapper for autodock_convert that takes a single structure folder path as input

    # parse filenames and filetypes from folder path
    folder_name = folder.split('/')
    folder_name = folder_name[len(folder_name) - 1]
    ligand_file = [filename for filename in os.listdir(folder) if 'ligand.pdb' in filename][0]
    ligand_filepath = folder + '/' + ligand_file
    receptor_file = [filename for filename in os.listdir(folder) if 'receptor.pdb' in filename or 'pocket.pdb' in filename][0]
    receptor_filepath = folder + '/' + receptor_file
    filetype = ligand_file.split('.')[1]

    # convert using autodock_convert function
    autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file, filetype)

def batch_convert_to_pdbqt(druglike_structure_folders, num_threads): # wrapper for multiprocessing convert_to_pdbqt function
    threads = min(num_threads, len(druglike_structure_folders))

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(convert_to_pdbqt, druglike_structure_folders), total=len(druglike_structure_folders)))

def move_problematic_structures(fatal_error_list, druglike_files_path, destination_path, problem_path): # copy structures in fatal_error_list to location for manual inspection
    fatal_error_list = list(set(list(fatal_error_list)))
    with tqdm(total=len(fatal_error_list)) as pbar:
        for folder in fatal_error_list:
            shutil.rmtree(f'{destination_path}{folder}')
            shutil.copytree(f'{druglike_files_path}{folder}', f'{problem_path}{folder}')
            pbar.update(1)

def parse_args(args): # parse CLI user inputs

    druglike_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-suc') + 1]

    problem_path = args[args.index('-prob') + 1]

    prep_ligand_command = args[args.index('-prep_lig') + 1]

    prep_protein_command = args[args.index('-prep_prot') + 1]

    num_threads = args[args.index('-threads') + 1]

    return druglike_files_path, destination_path, problem_path, prep_ligand_command, prep_protein_command, num_threads

def main(): # run script using CLI

    # temporary fix for global variable error
    global destination_path
    global prep_ligand_command
    global prep_protein_command

    druglike_files_path, destination_path, problem_path, prep_ligand_command, prep_protein_command, num_threads = parse_args(sys.argv)

    druglike_structure_folders = [(druglike_files_path + folder) for folder in os.listdir(druglike_files_path)]

    druglike_structure_folders = druglike_structure_folders

    batch_convert_to_pdbqt(druglike_structure_folders, num_threads)

    move_problematic_structures(fatal_error_list, druglike_files_path, destination_path, problem_path)

if __name__ == '__main__':
    main()
