import os
import oddt
import subprocess
from tqdm import tqdm

MOAD_fatal_error_list = list()

pdb_files_path = ''
moad_files_path = '/home/milesm/Dissertation/Data/Parsed/Non_redundant/Binding_MOAD/'
dude_files_path = ''

prep_ligand_command = '/home/milesm/Dissertation/Third_Party_Code/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
prep_protein_command = '/home/milesm/Dissertation/Third_Party_Code/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'

destination_path = '/home/milesm/Dissertation/Data/PDBQT/Non_redundant/'

example_ligand_file = '/home/milesm/Desktop/4yki_ligand.pdb'
example_protein_file = '/home/milesm/Desktop/7cpa/7cpa_protein.pdb'
example_destination_path = '/home/milesm/Desktop/7cpa/'

#pdb_structure_folders = [(pdb_files_path + folder) for folder in os.listdir(pdb_files_path)]
moad_structure_folders = [(moad_files_path + folder) for folder in os.listdir(moad_files_path)]
#dude_structure_folders = [(dude_files_path + folder) for folder in os.listdir(dude_files_path)]

def repair_stacked_ligand(ligand_filepath, ligand_file):
    ligand_raw_text = open(ligand_filepath, 'r').read()
    atoms = ligand_raw_text.split('\n')
    if atoms[0].strip().lstrip() == 'END':
        print(f'{ligand_file} is an empty file!')
        return None
    else:
        end_of_molecule = [idx for idx, atom in enumerate(atoms) if 'TER' in atom][0] + 1
        molecule = '\n'.join(atoms[:end_of_molecule])
        repaired_ligand_name = ligand_file.split('.')[0] + '_REPAIRED.pdb'
        ligand_base_path = ligand_filepath.split('/')
        ligand_base_path = ligand_base_path[:len(ligand_base_path) - 1]
        ligand_base_path = '/'.join(ligand_base_path)
        repaired_ligand_filepath = f'{ligand_base_path}/{repaired_ligand_name}'
        ligand_base_path_files = os.listdir(ligand_base_path)
        repair_check = [file for file in ligand_base_path_files if 'REPAIRED' in file]
        if len(repair_check) == 0:
            with open(repaired_ligand_filepath, 'a+') as repaired_ligand:
                repaired_ligand.write(molecule)
                repaired_ligand.close()
        else:
            print(f'Found pre-repaired ligand from a previous run: "{repair_check[0]}"')
        return repaired_ligand_filepath

def autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file):
    try:
        os.mkdir(f'{destination_path}{folder_name}')
    except FileExistsError:
        pass
    print(folder_name)
    print('Preparing ligand...')
    try:
        output = subprocess.check_output(f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_file}qt', shell=True, stderr=subprocess.STDOUT)
    except:
        try:
            print('REPAIRING LIGAND...')
            repaired_ligand_filepath = repair_stacked_ligand(ligand_filepath, ligand_file)
            if repaired_ligand_filepath is not None:
                output = subprocess.check_output(f'{prep_ligand_command} -l {repaired_ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_file}qt', shell=True)
                print('REPAIR SUCCESSFUL!')
        except:
            MOAD_fatal_error_list.append(folder_name)
            if os.path.isfile('MOAD_fatal_error_list.txt'):
                os.remove('MOAD_fatal_error_list.txt')
            with open('MOAD_fatal_error_list.txt','a+') as error_list:
                error_list.write(str(MOAD_fatal_error_list))
                error_list.close()
            print(f'FATAL PROBLEM WITH {folder_name}: Added to list and skipping...')

    print('Preparing protein...')
    output = subprocess.check_output(f'{prep_protein_command} -r {receptor_filepath} -A hydrogens -o {destination_path}{folder_name}/{receptor_file}qt -U waters', shell=True)


def batch_convert_to_pdbqt(structure_folders, database_name, destination_path):
    if database_name == 'moad':
        destination_path = destination_path + 'Binding_MOAD/'
        try:
            os.mkdir(destination_path)
        except FileExistsError:
            pass
        with tqdm(total=len(structure_folders)) as pbar:
            for folder in structure_folders:
                folder_name = folder.split('/')
                folder_name = folder_name[len(folder_name) - 1]
                files_in_folder = os.listdir(folder)
                ligand_file = [filename for filename in files_in_folder if 'ligand.pdb' in filename and '.pdbqt' not in filename][0]
                ligand_filepath = folder + '/' + ligand_file
                receptor_file = [filename for filename in files_in_folder if 'protein.pdb' in filename and '.pdbqt' not in filename][0]
                receptor_filepath = folder + '/' + receptor_file
                autodock_convert(folder_name, destination_path, prep_ligand_command, prep_protein_command, ligand_filepath, ligand_file, receptor_filepath, receptor_file)
                pbar.update(1)

batch_convert_to_pdbqt(moad_structure_folders, 'moad', destination_path)
