import os
import oddt
import subprocess
from tqdm import tqdm
import shutil
import sys

fatal_error_list = list()

def autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file, filetype):
    try:
        os.mkdir(f'{destination_path}{folder_name}')
        print(f'Made directory {destination_path}{folder_name}')
    except FileExistsError:
        pass
    print(folder_name)
    if filetype == 'mol2':
        ligand_outfile = ligand_file.replace('mol2','pdb')
    else:
        ligand_outfile = ligand_file
    print('Preparing ligand...')
    try:
        output = subprocess.check_output(f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_outfile}qt -U nphs', shell=True, stderr=subprocess.STDOUT)
    except:
        fatal_error_list.append(folder_name)
        if os.path.isfile('fatal_error_list.txt'):
            os.remove('fatal_error_list.txt')
        with open('fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list))
            error_list.close()
        print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

    print('Preparing protein...')
    try:
        output = subprocess.check_output(f'{prep_protein_command} -r {receptor_filepath} -A hydrogens -o {destination_path}{folder_name}/{receptor_file}qt -U nphs', shell=True)
    except:
        fatal_error_list.append(folder_name)
        if os.path.isfile('fatal_error_list.txt'):
            os.remove('fatal_error_list.txt')
        with open('fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list))
            error_list.close()
        print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

def batch_convert_to_pdbqt(structure_folders, destination_path, prep_ligand_command, prep_protein_command):
    try:
        os.mkdir(destination_path)
    except FileExistsError:
        pass
    with tqdm(total=len(structure_folders)) as pbar:
        for folder in structure_folders:
            print(folder)
            folder_name = folder.split('/')
            folder_name = folder_name[len(folder_name) - 1]
            ligand_file = [filename for filename in os.listdir(folder) if 'ligand.pdb' in filename][0]
            ligand_filepath = folder + '/' + ligand_file
            receptor_file = [filename for filename in os.listdir(folder) if 'receptor.pdb' in filename or 'pocket.pdb' in filename][0]
            receptor_filepath = folder + '/' + receptor_file
            filetype = ligand_file.split('.')[1]
            autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file, filetype)
            pbar.update(1)

def move_problematic_structures(fatal_error_list, druglike_files_path, destination_path, problem_path):
    fatal_error_list = list(set(list(fatal_error_list)))
    with tqdm(total=len(fatal_error_list)) as pbar:
        for folder in fatal_error_list:
            shutil.rmtree(f'{destination_path}{folder}')
            shutil.copytree(f'{druglike_files_path}{folder}', f'{problem_path}{folder}')
            pbar.update(1)

def parse_args(args):

    druglike_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-suc') + 1]

    problem_path = args[args.index('-prob') + 1]

    prep_ligand_command = args[args.index('-prep_lig') + 1]

    prep_protein_command = args[args.index('-prep_prot') + 1]

    return druglike_files_path, destination_path, problem_path, prep_ligand_command, prep_protein_command

def main():

    druglike_files_path, destination_path, problem_path, prep_ligand_command, prep_protein_command = parse_args(sys.argv)

    druglike_structure_folders = [(druglike_files_path + folder) for folder in os.listdir(druglike_files_path)]

    druglike_structure_folders = druglike_structure_folders

    batch_convert_to_pdbqt(druglike_structure_folders, destination_path, prep_ligand_command, prep_protein_command)

    move_problematic_structures(fatal_error_list, druglike_files_path, destination_path, problem_path)

if __name__ == '__main__':
    main()
