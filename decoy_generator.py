import oddt

example_active = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/actives/CHEMBL190.pdbqt'
destination_path = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/'

def convert_pdbqt_to_smiles(pdbqt_filepath, destination_path):
    pdbqt_filename = pdbqt_filepath.split('/')
    pdbqt_filename = pdbqt_filename[len(pdbqt_filename) - 1]
    pdb_code = pdbqt_filename.split('.')[0]
    pdbqt = next(oddt.toolkits.ob.readfile('pdbqt', pdbqt_filepath))
    pdbqt.write('smi',f'{destination_path}{pdb_code}.txt', overwrite=True)
    with open(f'{destination_path}{pdb_code}.txt', 'r') as smi_file:
        smile = smi_file.read().split(destination_path)[0]
    return (smile + ' ' + pdb_code)

def convert_pdb_to_smiles(pdb_filepath, destination_path):
    pdb_filename = pdb_filepath.split('/')
    pdb_filename = pdb_filename[len(pdb_filename) - 1]
    pdb_code = pdb_filename.split('.')[0]
    pdb = next(oddt.toolkits.ob.readfile('pdb', pdb_filepath))
    pdb.write('smi',f'{destination_path}{pdb_code}.txt', overwrite=True)
    with open(f'{destination_path}{pdb_code}.txt', 'r') as smi_file:
        smile = smi_file.read().split(destination_path)[0]
    return (smile + ' ' + pdb_code)

def convert_mol2_to_smiles(mol2_filepath, destination_path):
    mol2_filename = mol2_filepath.split('/')
    mol2_filename = mol2_filename[len(mol2_filename) - 1]
    pdb_code = mol2_filename.split('.')[0]
    mol2 = next(oddt.toolkits.ob.readfile('mol2', mol2_filepath))
    mol2.write('smi',f'{destination_path}{pdb_code}.txt', overwrite=True)
    with open(f'{destination_path}{pdb_code}.txt', 'r') as smi_file:
        smile = smi_file.read().split(destination_path)[0]
    return (smile + ' ' + pdb_code)

smile = convert_pdbqt_to_smiles(example_active, destination_path)
print(smile)
