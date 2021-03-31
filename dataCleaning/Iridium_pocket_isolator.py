'''
Loops through Iridium dataset pdb structures and isolates
a receptor.pdb file from the structure which includes any
residues with atoms within a user defined distance
(cutoff_thresh) of the crystal ligand
'''

import os
from Bio.PDB import *
from tqdm import tqdm
import pandas as pd
import oddt
import sys
from warnings import filterwarnings
from tqdm import tqdm
from itertools import chain as iterchain

class pocket_selector(Select): # selector class for BioPython which saves all residues in supplied 'pocket_residues' list

    def __init__(self, pocket_residues):
        self.pocket_residues = pocket_residues

    def accept_residue(self, residue): # accept residues if unique identifier present in 'pocket_residues' list
        return True if (str(residue.id[1]) + residue.get_parent().id) in self.pocket_residues else False

def define_pocket(structure, chain_residue, cutoff_thresh): # calculate cuboid coordinates around ligand with padding of user defined 'cutoff_thresh' angstroms

    # construct list of coordinates of all atoms in ligand
    atom_coord_data = list()
    for atom in chain_residue.get_atoms():
        atom_coord_data.append(atom.get_coord())

    # turn coordinates into dataframe
    atom_df = pd.DataFrame(atom_coord_data, columns=['X','Y','Z'])

    # calculate min and max of x, y, z coordinates and pad with 'cutoff_thresh' angstroms
    min_x, max_x = (min(atom_df.X) - cutoff_thresh), (max(atom_df.X) + cutoff_thresh)
    min_y, max_y = (min(atom_df.Y) - cutoff_thresh), (max(atom_df.Y) + cutoff_thresh)
    min_z, max_z = (min(atom_df.Z) - cutoff_thresh), (max(atom_df.Z) + cutoff_thresh)
    return (min_x, min_y, min_z, max_x, max_y, max_z)

def pocket_check(x, y, z, min_x, min_y, min_z, max_x, max_y, max_z): # check if an atom is within the cuboid produced by 'define_pocket' function
    if min_x < x < max_x and min_y < y < max_y and min_z < z < max_z:
        return True
    else:
        return False

def get_ligand_information(ligand): # construct dataframe of ligand residues id, number, insertion code, chain letter

    # make and return dataframe
    ligand_residue_data = list()
    chain_data = list()
    for residue in ligand.get_residues():
        ligand_residue_data.append(residue.id)
        chain_data.append(residue.get_parent().id)
    ligand_df = pd.DataFrame(ligand_residue_data, columns=['ID','SEQ','INS'])
    ligand_df['CHAIN'] = chain_data
    return ligand_df

def broken_ligand(ligand_filepath, file_format, ligand): # test for separate ligands stored as one chain

    ligand_information = get_ligand_information(ligand)

    # pass test if only single residue in ligand
    if len(ligand_information) == 1:
        return False
    else:

        # iterate through atoms and bonds to populate bonds_dictionary
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

        # check for continuous shared bonds between residues
        shared_bonds = None
        for index, item in enumerate(residue_identifiers):
            try:
                shared_bonds = len(list(set(bonds_dictionary[residue_identifiers[index]]).intersection(bonds_dictionary[residue_identifiers[index + 1]])))
            except IndexError:
                pass

        # fail the test if no shared bonds found between residues else pass
        if shared_bonds == 0:
            return True
        else:
            return False

def isolate_pocket(filename, success_destination_path, problem_destination_path, cutoff_thresh): # save all residues within cutoff_thresh of cuboid around ligand

    # default unless problem found
    keep_structure = True

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    # define variables from filename
    ligand_path = f'{filename}_deposited_refined_lig.sdf'
    mol = next(oddt.toolkits.ob.readfile('sdf', ligand_path))
    mol.write('pdb', f'{filename}_deposited_refined_lig.pdb', overwrite=True)
    ligand_path = f'{filename}_deposited_refined_lig.pdb'
    protein_path = f'{filename}_deposited_refined_prot.pdb'

    # define pdb code
    pdb_code = filename.split('/')[len(filename.split('/')) - 1]

    # load the pdb structure
    parser = PDBParser()
    protein = parser.get_structure(pdb_code, protein_path)
    ligand = parser.get_structure('ligand', ligand_path)

    # check for amino acids in ligand
    amino_acids_in_ligand = [residue.get_resname() for residue in ligand.get_residues() if residue.get_resname() in amino_acids]
    if len(amino_acids_in_ligand) != 0:
        print('Peptide Ligand detected!')
        keep_structure = False

    # check for breaks in ligand
    elif broken_ligand(ligand_path, 'pdb', ligand):
        print('Found broken ligand!')
        keep_structure = False

    else:

        # define dict for populating with cuboid dimensions for each residue in the ligand
        ligand_dimensions = dict()

        # define list for populating with residues found within user defined cutoff_thresh
        pocket_residues = list()

        # loop through residues in ligand and get surrounding cuboid dimensions
        for residue in ligand.get_residues():
            residue_unique_id = str(residue.id[1]) + residue.get_parent().id
            ligand_pocket_dimensions = define_pocket(ligand, residue, cutoff_thresh)
            ligand_dimensions[residue.id[0].replace('H_','')] = ligand_pocket_dimensions

        # loop through residues in protein and add to 'pocket_residues' if they have any atoms inside cuboid cutoff_thresh
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

    # save the receptor and ligand files if structure is being kept
    if keep_structure:
        try:
            os.mkdir(f'{success_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_receptor.pdb', pocket_selector(pocket_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_protein.pdb')
        mol.write('pdb', f'{success_destination_path}{pdb_code}/{pdb_code}_ligand.pdb', overwrite=True)
        print('Saved files!')

    # save copy of the pdb file if structure has problems
    else:
        try:
            os.mkdir(f'{problem_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        print('Bad Structure - Saving to problem directory and skipping..')
        io.save(f'{problem_destination_path}{pdb_code}/{pdb_code}.pdb')
        pass

def parse_args(args): # parse CLI user inputs

    structure_location = args[args.index('-loc') + 1]

    success_destination_path = args[args.index('-suc') + 1]

    problem_destination_path = args[args.index('-prob') + 1]

    cutoff_thresh = float(args[args.index('-cutoff') + 1])

    return structure_location, success_destination_path, problem_destination_path, cutoff_thresh

def main(): # run script using CLI

    structure_location, success_destination_path, problem_destination_path, cutoff_thresh = parse_args(sys.argv)

    structure_files = os.listdir(structure_location)

    structures = list(set(list([(structure_location + file.split('_')[0]) for file in structure_files])))

    # loop through structures and isolate pocket files
    with tqdm(total=len(structures)) as pbar:
        for structure_file in structures:
            isolate_pocket(structure_file, success_destination_path, problem_destination_path, cutoff_thresh)
            pbar.update(1)

if __name__ == '__main__':
    main()
