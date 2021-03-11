# this script should load up a full pdb or biounit file into biopython.
# Then get the coordinates of the ligand atoms, and define a box around the ligand
# Then check which atoms are inside the box by their coordinates, and get a list of the residues these atoms belong to
# Then it should save all the residues that are in that list in a receptor file, and save the ligand in a ligand file

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

class ligand_selector(Select):

    def __init__(self, keep_ligand_residues):
        self.keep_ligand_residues = keep_ligand_residues

    def accept_residue(self, residue):
        # select residues from identified ligand
        return True if (str(residue.id[1]) + residue.get_parent().id) in self.keep_ligand_residues else False

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

def get_ligand_information(structure, valid_ligand_ids, valid_chain_id):
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

def broken_ligand_check(filepath, file_format, valid_ligand_ids, valid_chain_id, structure):
    ligand_information = get_ligand_information(structure, valid_ligand_ids, valid_chain_id)
    if len(ligand_information) == 1:
        return True
    else:
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
        shared_bonds = None
        for index, item in enumerate(residue_identifiers):
            try:
                shared_bonds = len(list(set(bonds_dictionary[residue_identifiers[index]]).intersection(bonds_dictionary[residue_identifiers[index + 1]])))
            except IndexError:
                pass
        if shared_bonds == 0:
            return False
        else:
            return True

def isolate_pocket_and_ligand(filename, success_destination_path, problem_destination_path, datafile, cutoff_thresh, exclusion_thresh):

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

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
    valid_ligand_chain = [ligand.split(':')[1] for ligand in list(ligand_data.LigandID)][0]
    valid_ligand_ids = [ligand.split(':')[0].split(' ') for ligand in list(ligand_data.LigandID)]

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

    if not broken_ligand_check(filename, 'pdb', valid_ligand_ids, valid_ligand_chain, structure):
        keep_structure = False
        print('BROKEN LIGAND')
        pass

    else:

        ligand_dimensions = dict()

        exclusion_dimensions = dict()

        pocket_residues = list()

        keep_ligand_residues = list()

        ligand_count = 1

        model = structure[0]
        chain = model[valid_ligand_chain]

        for chain_residue in chain.get_residues():
            if chain_residue.id[0].replace('H_','') in valid_ligand_ids[0]:
                residue_unique_id = str(chain_residue.id[1]) + chain_residue.get_parent().id
                ligand_pocket_dimensions = define_pocket(structure, chain_residue, cutoff_thresh)
                ligand_exclusion = define_pocket(structure, chain_residue, exclusion_thresh)
                ligand_dimensions[chain_residue.id[0].replace('H_','')] = ligand_pocket_dimensions
                exclusion_dimensions[chain_residue.id[0].replace('H_','')] = ligand_exclusion
                keep_ligand_residues.append(residue_unique_id)

        if len(ligand_dimensions) == 0:
            # only true if no ligand contains no HETATMS - second check for peptide ligands stored as single letter amino acid codes
            print('Peptide Ligand detected!')
            keep_structure = False

        for residue in structure.get_residues():
            for atom in residue.get_atoms():
                x, y, z = atom.get_coord()
                for exclusion in exclusion_dimensions.values():
                    if pocket_check(x, y, z, *exclusion):
                        if 'H' in atom.get_parent().id[0]:
                            if atom.get_parent().id[0].replace('H_','') in valid_ligand_ids[0]:
                                if str(atom.get_parent().get_parent().id) == valid_ligand_chain:
                                    pass
                                elif atom.get_parent().id[0].replace('H_','') in protein_part_ids:
                                    pass
                                else:
                                    keep_structure = False


                for dimension in ligand_dimensions.values():
                    if pocket_check(x, y, z, *dimension):
                        residue_unique_id = str(atom.get_parent().id[1]) + atom.get_parent().get_parent().id
                        if atom.get_parent().id[0] == 'W':
                            pass
                        elif 'H' in atom.get_parent().id[0]:
                            if atom.get_parent().id[0].replace('H_','').strip().lstrip() in valid_ligand_ids[0]:
                                pass
                            elif atom.get_parent().id[0].replace('H_','').strip().lstrip() in protein_part_ids:
                                pocket_residues.append(residue_unique_id)
                        elif atom.get_parent().id[0] == ' ':
                            pocket_residues.append(residue_unique_id)

    if keep_structure:
        try:
            os.mkdir(f'{success_destination_path}{pdb_code}')
        except FileExistsError:
            pass
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_receptor.pdb', pocket_selector(pocket_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}_ligand.pdb', ligand_selector(keep_ligand_residues))
        io.save(f'{success_destination_path}{pdb_code}/{pdb_code}.pdb')
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

    binding_MOAD_datafile_path = args[args.index('-ref') + 1]

    cutoff_thresh = float(args[args.index('-cutoff') + 1])

    exclusion_thresh = float(args[args.index('-exclusion') + 1])

    return structure_location, success_destination_path, problem_destination_path, binding_MOAD_datafile_path, cutoff_thresh, exclusion_thresh

def main():

    structure_location, success_destination_path, problem_destination_path, binding_MOAD_datafile_path, cutoff_thresh, exclusion_thresh = parse_args(sys.argv)

    structure_folders = os.listdir(structure_location)

    structures = [(structure_location + folder) for folder in structure_folders]

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

    with tqdm(total=len(structures)) as pbar:
        for structure_file in structures:
            isolate_pocket_and_ligand(structure_file, success_destination_path, problem_destination_path, binding_df, cutoff_thresh, exclusion_thresh)
            pbar.update(1)

if __name__ == '__main__':
    main()
