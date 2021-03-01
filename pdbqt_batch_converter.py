import os
import oddt

pdb_files_path = ''
moad_files_path = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/'
dude_files_path = ''

prep_ligand_command = 'autodocktools_prep_ligand'
prep_protein_command = 'autodocktools_prep_receptor'

destination_path = '/home/milesm/Dissertation/Data/PDBQT/'

example_ligand_file = '/home/milesm/Desktop/7cpa/7cpa_ligand.pdb'
example_protein_file = '/home/milesm/Desktop/7cpa/7cpa_protein.pdb'
example_destination_path = '/home/milesm/Desktop/7cpa/'

#pdb_structure_folders = [(pdb_files_path + folder) for folder in os.listdir(pdb_files_path)]
moad_structure_folders = [(moad_files_path + folder) for folder in os.listdir(moad_files_path)]
#dude_structure_folders = [(dude_files_path + folder) for folder in os.listdir(dude_files_path)]

def repair_stacked_ligand(ligand_filepath):
    ligand_raw_text = open(ligand_filepath, 'r')
    print(ligand_raw_text.split('TER'))

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
                print(ligand_file)
                print(receptor_file)
                try:
                    os.mkdir(f'{destination_path}{folder_name}')
                except FileExistsError:
                    pass
                sp = subprocess.Popen(["/bin/bash", "-i", "-c", f'{prep_ligand_command} -l {ligand_filepath} -A hydrogens -o {destination_path}Binding_MOAD/{folder_name}/{ligand_file}qt'])
                ligand_output = sp.communicate()[0]
                print(ligand_output)
                #if 'AttributeError' in ligand_output:
                #    repair_stacked_ligand(ligand_filepath)
                sp = subprocess.Popen(["/bin/bash", "-i", "-c", f'{prep_protein_command} -r {ligand_filepath} -A hydrogens -o {destination_path}Binding_MOAD/{folder_name}/{receptor_file}qt -U waters'])
                receptor_output = sp.communicate()[0]
                input('')


def convert_to_pdbqt(pdb_filepath, destination_path):
    pdb_filename = pdb_filepath.split('/')
    pdb_filename = pdb_filename[len(pdb_filename) - 1]
    pdb = next(oddt.toolkits.ob.readfile('pdb', pdb_filepath))
    pdb.addh()
    pdb.calccharges()
    pdb.write('pdbqt',f'{destination_path}{pdb_filename}qt', overwrite=True)



convert_to_pdbqt(example_ligand_file, example_destination_path)
convert_to_pdbqt(example_protein_file, example_destination_path)
