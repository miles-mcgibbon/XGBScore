import os
from Bio.PDB import *
from tqdm import tqdm
import pandas as pd
import oddt
import sys
from warnings import filterwarnings
from tqdm import tqdm
from itertools import chain as iterchain

# ignore warnings about incomplete PDB structures, as Binding MOAD supplies truncated biounits
filterwarnings('ignore')

class pocket_selector(Select):

    def __init__(self, pocket_residues):
        self.pocket_residues = pocket_residues

    def accept_residue(self, residue):
        # select residues from identified ligand
        return True if (str(residue.id[1]) + residue.get_parent().id) in self.pocket_residues else False

def define_pocket(structure, chain_residue, boundary):
    atom_coord_data = list()
    for atom in chain_residue.get_atoms():
        atom_coord_data.append(atom.get_coord())
    atom_df = pd.DataFrame(atom_coord_data, columns=['X','Y','Z'])
    min_x, max_x = (min(atom_df.X) - boundary), (max(atom_df.X) + boundary)
    min_y, max_y = (min(atom_df.Y) - boundary), (max(atom_df.Y) + boundary)
    min_z, max_z = (min(atom_df.Z) - boundary), (max(atom_df.Z) + boundary)
    return (min_x, min_y, min_z, max_x, max_y, max_z)

def pocket_check(x, y, z, min_x, min_y, min_z, max_x, max_y, max_z):
    if min_x < x < max_x and min_y < y < max_y and min_z < z < max_z:
        return True
    else:
        return False

def get_ligand_information(ligand):
    ligand_residue_data = list()
    chain_data = list()
    for residue in ligand.get_residues():
        ligand_residue_data.append(residue.id)
        chain_data.append(residue.get_parent().id)
    ligand_df = pd.DataFrame(ligand_residue_data, columns=['ID','SEQ','INS'])
    ligand_df['CHAIN'] = chain_data
    return ligand_df

def broken_ligand(ligand_filepath, file_format, ligand):
    ligand_information = get_ligand_information(ligand)
    if len(ligand_information) == 1:
        return False
    else:
        mol = next(oddt.toolkits.ob.readfile(file_format, ligand_filepath))
        bonds_dictionary = dict()
        residue_identifiers = list()
        for residue in mol.residues:
            residue_identifier = str(residue.idx)
            residue_identifiers.append(residue_identifier)
            atoms_in_residue = list()
            for atom in residue.atoms:
                for bond in atom.bonds:
                    for atom in bond.atoms:
                        atoms_in_residue.append(atom.idx)
            bonds_dictionary[residue_identifier] = atoms_in_residue
        shared_bonds = None
        for index, item in enumerate(residue_identifiers):
            try:
                shared_bonds = len(list(set(bonds_dictionary[residue_identifiers[index]]).intersection(bonds_dictionary[residue_identifiers[index + 1]])))
            except IndexError:
                pass
        if shared_bonds == 0:
            return True
        else:
            return False

def isolate_pocket_and_ligand(filename, success_destination_path, problem_destination_path, cutoff_thresh):

    keep_structure = True

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    ligand = [f'{filename}/{file}' for file in os.listdir(filename) if 'ligand.mol2' in file][0]
    ext = ligand.split('.')[1]
    print(ligand)
    print(ext)
    mol = next(oddt.toolkits.ob.readfile(ext, ligand))
    mol.write('pdb', f"{ligand.split('.')[0]}.pdb", overwrite=True)
    ligand_path = [f'{filename}/{file}' for file in os.listdir(filename) if 'ligand.pdb' in file][0]
    protein_path = [f'{filename}/{file}' for file in os.listdir(filename) if 'protein' in file][0]

    # define pdb code
    pdb_code = filename.split('.')[0].split('/')[filename.count('/')]
    ext = ligand_path.split('.')[1]

    # load the pdb structure
    parser = PDBParser()
    protein = parser.get_structure(pdb_code, protein_path)
    ligand = parser.get_structure('ligand', ligand_path)

    amino_acids_in_ligand = [residue.get_resname() for residue in ligand.get_residues() if residue.get_resname() in amino_acids]
    if len(amino_acids_in_ligand) != 0:
        print('Peptide Ligand detected!')
        keep_structure = False

    elif broken_ligand(ligand_path, 'pdb', ligand):
        print('Found broken ligand!')
        keep_structure = False

    else:

        ligand_dimensions = dict()

        pocket_residues = list()

        for residue in ligand.get_residues():
            residue_unique_id = str(residue.id[1]) + residue.get_parent().id
            ligand_pocket_dimensions = define_pocket(ligand, residue, cutoff_thresh)
            ligand_dimensions[residue.id[0].replace('H_','')] = ligand_pocket_dimensions


        for residue in protein.get_residues():
            for atom in residue.get_atoms():
                x, y, z = atom.get_coord()
                for dimension in ligand_dimensions.values():
                    if pocket_check(x, y, z, *dimension):
                        if residue.get_resname() != 'HOH':
                            residue_unique_id = str(atom.get_parent().id[1]) + atom.get_parent().get_parent().id
                            pocket_residues.append(residue_unique_id)

    # set the structure for saving
    io = PDBIO()
    io.set_structure(protein)

    if keep_structure:
        try:
            os.mkdir(f'{success_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_receptor.pdb', pocket_selector(pocket_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_protein.pdb')
        mol.write('pdb', f'{success_destination_path}{pdb_code}/{pdb_code}_ligand.pdb', overwrite=True)
        print('Saved files!')
    else:
        try:
            os.mkdir(f'{problem_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        print('Bad Structure - Saving to problem directory and skipping..')
        io.save(f'{problem_destination_path}{pdb_code}/{pdb_code}.pdb')
        pass

def parse_args(args):

    structure_location = args[args.index('-loc') + 1]

    success_destination_path = args[args.index('-suc') + 1]

    problem_destination_path = args[args.index('-prob') + 1]

    cutoff_thresh = float(args[args.index('-cutoff') + 1])

    return structure_location, success_destination_path, problem_destination_path, cutoff_thresh

def main():

    structure_location, success_destination_path, problem_destination_path, cutoff_thresh = parse_args(sys.argv)

    structure_folders = os.listdir(structure_location)

    structures = [(structure_location + folder) for folder in structure_folders]

    with tqdm(total=len(structures)) as pbar:
        for structure_file in structures:
            if 'index' in structure_file or 'readme' in structure_file:
                pass
            else:
                isolate_pocket_and_ligand(structure_file, success_destination_path, problem_destination_path, cutoff_thresh)
                pbar.update(1)

if __name__ == '__main__':
    main()
