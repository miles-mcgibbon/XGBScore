import os
from Bio.PDB import *
from tqdm import tqdm
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

# load the pdb structure
parser = PDBParser()
structure = parser.get_structure("Test","/home/milesm/Desktop/1a0q.pdb")


# find all the het residues/non amino acids
residues = list(structure.get_residues())
ligands = [res for res in residues if 'H_' in res.id[0]]

# get the number of atoms all the ligands and store them in a dictionary
ligand_atoms = {}
for lig in ligands:
    ligand_atoms[len(lig.child_dict.values())] = lig.id


# get the longest ligand
active_ligand_id = ligand_atoms[max(ligand_atoms.keys())][0]

# save the ligand
io = PDBIO()
io.set_structure(structure)
for residue in structure.get_residues():
    io.save("ligand.pdb", ligand_selector(residue, active_ligand_id))

# save the protein
io.save("protein.pdb", protein_selector())
