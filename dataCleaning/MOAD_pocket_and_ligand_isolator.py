'''
Loops through Binding MOAD original pdb structures and
splits them into respective ligand.pdb and receptor.pdb
files; receptor.pdb file contains any residues with atoms
within a user defined distance (cutoff_thresh) of the
crystal ligand
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

class ligand_selector(Select): # selector class for BioPython which saves all residues in supplied 'keep_ligand_residues' list

    def __init__(self, keep_ligand_residues):
        self.keep_ligand_residues = keep_ligand_residues

    def accept_residue(self, residue): # accept residues if unique identifier present in 'keep_ligand_residues' list
        return True if (str(residue.id[1]) + residue.get_parent().id) in self.keep_ligand_residues else False

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

def get_ligand_information(structure, valid_ligand_ids, valid_chain_id):  # construct dataframe of ligand residues id, number, insertion code, chain letter
    ligand_residue_data = list()
    chain_data = list()
    for residue in structure.get_residues():
        chain_id = residue.get_parent().id
        if residue.id[0].replace('H_','') in valid_ligand_ids[0]:
            ligand_residue_data.append(residue.id)
            chain_data.append(chain_id)
    ligand_df = pd.DataFrame(ligand_residue_data, columns=['ID','SEQ','INS'])
    ligand_df['CHAIN'] = chain_data
    return ligand_df.loc[ligand_df.CHAIN == valid_chain_id]

def broken_ligand_check(filepath, file_format, valid_ligand_ids, valid_chain_id, structure): # test for separate ligands stored as one chain

    # pass test if only single residue in ligand
    ligand_information = get_ligand_information(structure, valid_ligand_ids, valid_chain_id)
    if len(ligand_information) == 1:
        return True
    else:

        # iterate through atoms and bonds to populate bonds_dictionary
        mol = next(oddt.toolkits.ob.readfile(file_format, filepath))
        bonds_dictionary = dict()
        residue_identifiers = list()
        for residue in mol.residues:
            if residue.name in valid_ligand_ids[0] and residue.chain.upper() == valid_chain_id:
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
            return False
        else:
            return True

def isolate_pocket_and_ligand(filename, success_destination_path, problem_destination_path, datafile, cutoff_thresh, exclusion_thresh): # save separate receptor and ligand pdb files

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    # default keep structure unless problem found
    keep_structure = True

    # define pdb code
    pdb_code = filename.split('.')[0].split('/')[filename.count('/')]
    ext = filename.split('.')[1].replace('bio','')

    # load the pdb structure
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, filename)

    # get the id of the active ligand
    pdb_code_data = datafile.loc[datafile.PDBCode.str.upper() == pdb_code.upper()]
    ligand_data = pdb_code_data.loc[pdb_code_data.Validity == 'valid']

    # get the chain id of one active ligand if multiple present
    valid_ligand_chain = [ligand.split(':')[1] for ligand in list(ligand_data.LigandID)][0]
    valid_ligand_ids = [ligand.split(':')[0].split(' ') for ligand in list(ligand_data.LigandID)]

    # discard structure if amino acids present in ligand
    if len(list(set(valid_ligand_ids[0]).intersection(amino_acids))) != 0:
        print('Peptide Ligand detected!')
        keep_structure = False

    # get the ids of any hetatms that are part of the protein
    pdb_code_data = datafile.loc[datafile.PDBCode.str.upper() == pdb_code.upper()]
    protein_parts = pdb_code_data.loc[pdb_code_data.Validity == 'Part of Protein']
    protein_part_ids = [part.split(':')[0].split(' ') for part in list(protein_parts.LigandID)]
    protein_part_ids = list(iterchain.from_iterable(protein_part_ids))
    if len(protein_part_ids) == 0:
        protein_part_ids.append(['NO KEY PROTEIN HETATMS'])



    # set the structure for saving
    io = PDBIO()
    io.set_structure(structure)

    # check ligand for residues with no shared bonds
    if not broken_ligand_check(filename, 'pdb', valid_ligand_ids, valid_ligand_chain, structure):
        keep_structure = False
        print('BROKEN LIGAND')
        pass

    else:

        # set up lists and dicts for populating
        ligand_dimensions = dict()
        exclusion_dimensions = dict()
        pocket_residues = list()
        keep_ligand_residues = list()

        model = structure[0]
        chain = model[valid_ligand_chain]

        # loop through residues in chosen active ligand chain
        for chain_residue in chain.get_residues():
            if chain_residue.id[0].replace('H_','') in valid_ligand_ids[0]:

                # if residue is in valid ligands make unique ID
                residue_unique_id = str(chain_residue.id[1]) + chain_residue.get_parent().id

                # calculate user defined receptor cuboid around the ligand
                ligand_pocket_dimensions = define_pocket(structure, chain_residue, cutoff_thresh)
                ligand_dimensions[chain_residue.id[0].replace('H_','')] = ligand_pocket_dimensions

                # calculate user defined exclusion cuboid around the ligand
                ligand_exclusion = define_pocket(structure, chain_residue, exclusion_thresh)
                exclusion_dimensions[chain_residue.id[0].replace('H_','')] = ligand_exclusion

                # add ligand residues to list for ligand selector class
                keep_ligand_residues.append(residue_unique_id)

        # second check for peptide ligands stored as single letter amino acid codes
        if len(ligand_dimensions) == 0:
            print('Peptide Ligand detected!')
            keep_structure = False

        # loop through protein residues
        for residue in structure.get_residues():
            for atom in residue.get_atoms():
                x, y, z = atom.get_coord()
                # check if any hetatms from a different chain/ligand are within the exclusion zone
                for exclusion in exclusion_dimensions.values():
                    if pocket_check(x, y, z, *exclusion):
                        if 'H' in atom.get_parent().id[0]:
                            if atom.get_parent().id[0].replace('H_','') in valid_ligand_ids[0]:
                                # ignore if hetatms are from the ligand itself
                                if str(atom.get_parent().get_parent().id) == valid_ligand_chain:
                                    pass
                                # ignore if hetatms are from bound metal ions or cofactors listed as part of protein
                                elif atom.get_parent().id[0].replace('H_','') in protein_part_ids:
                                    pass
                                # discard structure if hetatms found in exclusion zone that are not either of the above exceptions
                                else:
                                    keep_structure = False

                # find and record any protein residues with atoms inside the receptor cutoff
                for dimension in ligand_dimensions.values():
                    if pocket_check(x, y, z, *dimension):
                        residue_unique_id = str(atom.get_parent().id[1]) + atom.get_parent().get_parent().id
                        # ignore and remove water molecules
                        if atom.get_parent().id[0] == 'W':
                            pass
                        elif 'H' in atom.get_parent().id[0]:
                            # ignore ligand atoms
                            if atom.get_parent().id[0].replace('H_','').strip().lstrip() in valid_ligand_ids[0]:
                                pass
                            # save hetatms if listed as part of protein
                            elif atom.get_parent().id[0].replace('H_','').strip().lstrip() in protein_part_ids:
                                pocket_residues.append(residue_unique_id)
                        # save amino acid residues
                        elif atom.get_parent().id[0] == ' ':
                            pocket_residues.append(residue_unique_id)

    # save ligand and receptor pdb files if no problems found
    if keep_structure:
        try:
            os.mkdir(f'{success_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_receptor.pdb', pocket_selector(pocket_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_ligand.pdb', ligand_selector(keep_ligand_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}.pdb')
        print('Saved files!')

    # otherwise make copy of original pdb structure for manual verification
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

    binding_MOAD_datafile_path = args[args.index('-ref') + 1]

    cutoff_thresh = float(args[args.index('-cutoff') + 1])

    exclusion_thresh = float(args[args.index('-exclusion') + 1])

    return structure_location, success_destination_path, problem_destination_path, binding_MOAD_datafile_path, cutoff_thresh, exclusion_thresh

def main(): # run script using CLI

    structure_location, success_destination_path, problem_destination_path, binding_MOAD_datafile_path, cutoff_thresh, exclusion_thresh = parse_args(sys.argv)

    structure_folders = os.listdir(structure_location)

    structures = [(structure_location + folder) for folder in structure_folders]

    # load binding MOAD reference dataframe with forced headers and clean up formatting
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
    binding_df = pd.read_csv(binding_MOAD_datafile_path, names=force_columns)
    binding_df['PDBCode'] = binding_df['PDBCode'].ffill()

    # loop through structures and isolate receptor and ligand files
    with tqdm(total=len(structures)) as pbar:
        for structure_file in structures:
            isolate_pocket_and_ligand(structure_file, success_destination_path, problem_destination_path, binding_df, cutoff_thresh, exclusion_thresh)
            pbar.update(1)

if __name__ == '__main__':
    main()
