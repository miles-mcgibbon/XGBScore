import os
import oddt
import subprocess

pdb_files_path = ''
moad_files_path = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/'
dude_files_path = ''

prep_ligand_command = 'autodocktools_prep_ligand'
prep_protein_command = 'autodocktools_prep_receptor'

destination_path = '/home/milesm/Dissertation/Data/PDBQT/'

example_ligand_file = '/home/milesm/Desktop/4yki_ligand.pdb'
example_protein_file = '/home/milesm/Desktop/7cpa/7cpa_protein.pdb'
example_destination_path = '/home/milesm/Desktop/7cpa/'

#pdb_structure_folders = [(pdb_files_path + folder) for folder in os.listdir(pdb_files_path)]
moad_structure_folders = [(moad_files_path + folder) for folder in os.listdir(moad_files_path)]
#dude_structure_folders = [(dude_files_path + folder) for folder in os.listdir(dude_files_path)]

def repair_stacked_ligand(ligand_filepath):
    ligand_raw_text = open(ligand_filepath, 'r').read()
    atoms = ligand_raw_text.split('\n')
    end_of_molecule = [idx for idx, atom in enumerate(atoms) if 'TER' in atom][0] + 1
    molecule = '\n'.join(atoms[:end_of_molecule])
    os.remove(ligand_filepath)
    with open(ligand_filepath, 'a+') as repaired_ligand:
        repaired_ligand.write(molecule)
        repaired_ligand.close()


def batch_convert_to_pdbqt(structure_folders, database_name, destination_path):
    if database_name == 'moad':
        destination_path = destination_path + 'Binding_MOAD/'
        try:
            os.mkdir(destination_path)
        except FileExistsError:
            pass
        for folder in structure_folders:
            folder_name = folder.split('/')
            folder_name = folder_name[len(folder_name) - 1]
            if 'LIGAND_ERROR' in folder_name:
                pass
            else:
                files_in_folder = os.listdir(folder)
                ligand_file = [filename for filename in files_in_folder if 'ligand.pdb' in filename and '.pdbqt' not in filename][0]
                ligand_filepath = folder + '/' + ligand_file
                receptor_file = [filename for filename in files_in_folder if 'protein.pdb' in filename and '.pdbqt' not in filename][0]
                receptor_filepath = folder + '/' + receptor_file
                try:
                    os.mkdir(f'{destination_path}{folder_name}')
                except FileExistsError:
                    pass
                print('Preparing ligand...')
                sp = subprocess.Popen(["/bin/bash", "-i", "-c", f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_file}qt'])
                sp.communicate()[0]
                if os.path.isfile(f'{destination_path}{folder_name}/{ligand_file}qt'):
                    pass
                else:
                    print('Error - attempting repair...')
                    repair_stacked_ligand(ligand_filepath)
                    print('Preparing ligand...')
                    sp = subprocess.Popen(["/bin/bash", "-i", "-c", f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{ligand_file}qt'])
                    sp.communicate()[0]
                    if os.path.isfile(f'{destination_path}{folder_name}/{ligand_file}qt'):
                        pass
                    else:
                        print('Ligand error - could not resolve')
                print('Preparing protein...')
                sp = subprocess.Popen(["/bin/bash", "-i", "-c", f'{prep_protein_command} -r {ligand_filepath} -A hydrogens -o {destination_path}{folder_name}/{receptor_file}qt -U waters'])
                sp.communicate()[0]
                if os.path.isfile(f'{destination_path}{folder_name}/{receptor_file}qt'):
                    pass
                else:
                    print('Error')
                    # deal with exception
                print('Finished instance')
                input('')

batch_convert_to_pdbqt(moad_structure_folders, 'moad', destination_path)
