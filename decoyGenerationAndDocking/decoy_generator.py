import oddt
import os
from tqdm import tqdm

pdbqt_druglike_files_path = '/home/milesm/Dissertation/Data/PDBQT/Non_redundant_druglike/'

pdbqt_folders = [(pdbqt_druglike_files_path + folder) for folder in os.listdir(pdbqt_druglike_files_path)]

smiles_file = '/home/milesm/Dissertation/Data/PDBQT/smiles_for_decoys.txt'

def convert_file_to_smiles(filepath, filetype):
    filename = filepath.split('/')
    filename = filename[len(filename) - 1]
    pdb_code = filename.split('.')[0]
    molecule = next(oddt.toolkits.ob.readfile(filetype, filepath))
    return molecule.smiles + '    ' + pdb_code

with tqdm(total=len(pdbqt_folders)) as pbar:
    for filepath in pdbqt_folders:
        foldername = filepath.split('/')
        foldername = foldername[len(foldername) - 1]
        ligand_file = [filename for filename in os.listdir(filepath) if 'ligand.mol2' in filename or 'ligand.pdb' in filename][0]
        filetype = ligand_file.split('.')[1]
        ligand_filepath = filepath + '/' + ligand_file
        smile_line_to_write = convert_file_to_smiles(ligand_filepath, filetype)
        with open(smiles_file, 'a') as smiles_batch:
            smiles_batch.write(f'{smile_line_to_write}\n')
            smiles_batch.close()
        pbar.update(1)
