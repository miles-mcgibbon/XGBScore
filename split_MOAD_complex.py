import os
from Bio.PDB import *
from tqdm import tqdm
import pandas as pd
from warnings import filterwarnings

# ignore warnings about incomplete PDB structures, as Binding MOAD supplies truncated biounits
filterwarnings('ignore')

class ligand_selector(Select):

    def __init__(self, id_codes):
        self.id_codes = id_codes

    def accept_residue(self, residue):
        # select residues from identified ligand
        return True if residue.id[0].replace('H_','') in self.id_codes else False

class protein_selector(Select):

    def accept_residue(self, residue):
        # select protein residues
        return True if residue.id[0] == " " else False

def split_structure(filename, destination_path, datafile, make_dir, ignore_peptide):

    global ligand_errors

    # define pdb code
    pdb_code = filename.split('.')[0].split('/')[filename.count('/')]
    ext = filename.split('.')[1].replace('bio','')

    # load the pdb structure
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, filename)

    pdb_code_data = datafile.loc[datafile.PDBCode.str.upper() == pdb_code.upper()]
    ligand_data = pdb_code_data.loc[pdb_code_data.Validity == 'valid']
    valid_ligand_ids = [ligand.split(':')[0].split(' ') for ligand in list(ligand_data.LigandID)]

    # set the structure for saving
    io = PDBIO()
    io.set_structure(structure)

    if make_dir:
        try:
            os.mkdir(f'{destination_path}{pdb_code}')
        except FileExistsError:
            new_pdb_code = pdb_code + f'_{ext}'
            os.mkdir(f'{destination_path}{new_pdb_code}')
            pass

    for id_codes in valid_ligand_ids:
        io.save(f"{destination_path}{pdb_code}/{pdb_code}_ligand.pdb", ligand_selector(id_codes))

    io.save(f"{destination_path}{pdb_code}/{pdb_code}_protein.pdb", protein_selector())

complex_path = '/home/milesm/Dissertation/Data/Raw/Binding_MOAD/Extracted/BindingMOAD_2020/'
destination_path = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/'
binding_MOAD_datafile_path = '/home/milesm/Dissertation/Data/Raw/Binding_MOAD/Compressed/nr_bind.csv'

force_columns = ['Decimal',
                 'EntryString',
                 'PDBCode',
                 'LigandID',
                 'Validity',
                 'BindingDataType',
                 'Equals',
                 'BindingValue',
                 'BindingUnits',
                 'FormulaString',
                 'FormulaStringSpill']

complexes = os.listdir(complex_path)

ligand_errors = 0

binding_df = pd.read_csv(binding_MOAD_datafile_path, names=force_columns)

binding_df['PDBCode'] = binding_df['PDBCode'].ffill()

with tqdm(total=len(complexes)) as pbar:
    for complex in complexes:
        filename = complex_path + complex
        split_structure(filename, destination_path, binding_df, True, True)
        pbar.update(1)

print(f'Splitting completed with {ligand_errors} recorded ligand errors')
