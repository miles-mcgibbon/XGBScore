import os
from Bio.PDB import *
from tqdm import tqdm
import pandas as pd
from warnings import filterwarnings

# ignore warnings about incomplete PDB structures, as Binding MOAD supplies truncated biounits
filterwarnings('ignore')

class ligand_selector(Select):

    def __init__(self, residue, id_code):
        self.residue = residue
        self.id_code = id_code

    def accept_residue(self, residue):
        # select residues from identified ligand
        return True if residue.id[0] == self.id_code else False

class protein_selector(Select):

    def accept_residue(self, residue):
        # select protein residues
        return True if residue.id[0] == " " else False

def split_structure(filename, destination_path, datafile, make_dir, ignore_peptide):

    global ligand_errors

    # list of hetatm ids to ignore
    ignore_list = ['SO4','HOH']

    # define pdb code
    pdb_code = filename.split('.')[0].split('/')[filename.count('/')]

    # load the pdb structure
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, filename)


    # find all the het residues/non amino acids
    residues = list(structure.get_residues())
    ligands = [res for res in residues if 'H_' in res.id[0]]

    # get the number of atoms all the ligands and store them in a dictionary
    ligand_atoms = {}
    for lig in ligands:
        ligand_atoms[len(lig.child_dict.values())] = lig.id

    try:
        # get the longest ligands, excluding water and metal ions
        active_ligand_ids = [ligand_atoms[key][0] for key in ligand_atoms.keys() if key > 3 and ligand_atoms[key][0].replace('H_','') not in ignore_list]
        print(active_ligand_ids)

        # set the structure for saving
        io = PDBIO()
        io.set_structure(structure)

        if make_dir:
            try:
                os.mkdir(f'{destination_path}{pdb_code}')
            except FileExistsError:
                pass

            # save the ligands
            for residue in structure.get_residues():
                if residue.id[0] in active_ligand_ids:
                    io.save(f"{destination_path}{pdb_code}/{pdb_code}_ligand_{residue.id[0]}.pdb", ligand_selector(residue, residue.id[0]))

        # save the protein
        io.save(f"{destination_path}{pdb_code}/{pdb_code}_protein.pdb", protein_selector())
        input('')

    except ValueError:

        if ignore_peptide:
            pass
        else:
            # add flag to directory if no ligand found

            # set the structure for saving
            io = PDBIO()
            io.set_structure(structure)

            if make_dir:
                try:
                    os.mkdir(f'{destination_path}{pdb_code}_LIGAND_ERROR')
                except FileExistsError:
                    pass

            # save the protein
            io.save(f"{destination_path}{pdb_code}_LIGAND_ERROR/{pdb_code}_protein.pdb", protein_selector())

        ligand_errors = ligand_errors + 1

complex_path = '/home/milesm/Dissertation/Data/Raw/Binding_MOAD/Extracted/BindingMOAD_2020/'
destination_path = '/home/milesm/Dissertation/Data/Parsed/Binding_MOAD/'
binding_MOAD_datafile_path = '/home/milesm/Dissertation/Data/Raw/Binding_MOAD/Compressed/nr_bind.csv'

complexes = os.listdir(complex_path)

ligand_errors = 0

binding_df = pd.read_csv(binding_MOAD_datafile_path)

with tqdm(total=len(complexes)) as pbar:
    for complex in complexes:
        filename = complex_path + complex
        split_structure(filename, destination_path, binding_df, True, True)
        pbar.update(1)

print(f'Splitting completed with {ligand_errors} recorded ligand errors')
