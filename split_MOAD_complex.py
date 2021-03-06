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

    def __init__(self, id_codes):
        self.id_codes = id_codes

    def accept_residue(self, residue):
        # select protein residues
        if residue.id[0] == " ":
            return True
        elif residue.id[0] in self.id_codes:
            return True
        elif residue.id[0].replace('H_','') in self.id_codes:
            return True
        else:
            return False

def split_structure(filename, destination_path, datafile, make_dir, filter_peptides):

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    global ligand_errors
    global peptide_ligands

    # define pdb code
    pdb_code = filename.split('.')[0].split('/')[filename.count('/')]
    ext = filename.split('.')[1].replace('bio','')

    # load the pdb structure
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, filename)

    # get the id of the active ligand
    pdb_code_data = datafile.loc[datafile.PDBCode.str.upper() == pdb_code.upper()]
    ligand_data = pdb_code_data.loc[pdb_code_data.Validity == 'valid']
    valid_ligand_ids = [ligand.split(':')[0].split(' ') for ligand in list(ligand_data.LigandID)]

    # get the ids of any hetatms that are part of the protein
    pdb_code_data = datafile.loc[datafile.PDBCode.str.upper() == pdb_code.upper()]
    protein_parts = pdb_code_data.loc[pdb_code_data.Validity == 'Part of Protein']
    protein_part_ids = [part.split(':')[0].split(' ') for part in list(protein_parts.LigandID)]

    # set the structure for saving
    io = PDBIO()
    io.set_structure(structure)

    pdb_directory_name = pdb_code + '_' + ext

    peptide = False

    if filter_peptides:
        # check the ligand for amino acids
        for id_codes in valid_ligand_ids:
            for id_code in id_codes:
                id_code = id_code.strip().lstrip()
                if id_code in amino_acids:
                    peptide = True

    if peptide:
        peptide_ligands = peptide_ligands + 1
        pass
    else:
        if make_dir:
            try:
                os.mkdir(f'{destination_path}{pdb_directory_name}')
            except FileExistsError:
                pass
        for id_codes in valid_ligand_ids:
            io.save(f"{destination_path}{pdb_directory_name}/{pdb_directory_name}_ligand.pdb", ligand_selector(id_codes))
        if len(protein_part_ids) == 0:
            io.save(f"{destination_path}{pdb_directory_name}/{pdb_directory_name}_protein.pdb", protein_selector(id_codes))
        else:
            for id_codes in protein_part_ids:
                io.save(f"{destination_path}{pdb_directory_name}/{pdb_directory_name}_protein.pdb", protein_selector(id_codes))

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

peptide_ligands = 0

binding_df = pd.read_csv(binding_MOAD_datafile_path, names=force_columns)

binding_df['PDBCode'] = binding_df['PDBCode'].ffill()

with tqdm(total=len(complexes)) as pbar:
    for complex in complexes:
        filename = complex_path + complex
        split_structure(filename, destination_path, binding_df, True, True)
        pbar.update(1)

print(f'Splitting completed with {ligand_errors} recorded ligand errors')
print(f'Splitting completed with {peptide_ligands} recorded peptides that were skipped')
