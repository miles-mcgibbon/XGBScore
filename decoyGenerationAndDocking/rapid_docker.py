import os
from biopandas.mol2 import PandasMol2

example_crystal_ligand = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/crystal_ligand.mol2'
example_active = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/actives/CHEMBL190.pdbqt'
example_receptor = '/home/milesm/Dissertation/Data/Parsed/DUD-E/Separated/aa2ar/receptor.pdbqt'
docker_command = 'cd /home/milesm/Dissertation/Third_Party_Code/gwovina-1.0/build/linux/release && ./gwovina'

def get_coordinates(file, size):
    pmol = PandasMol2().read_mol2(file)
    x_min, x_max = min(pmol.df.x), max(pmol.df.x)
    x_center = x_min + abs(x_min-x_max)/2
    x_range = (abs(x_min-x_max)*size)
    y_min, y_max = min(pmol.df.y), max(pmol.df.y)
    y_center = y_min + abs(y_min-y_max)/2
    y_range = (abs(y_min-y_max)*size)
    z_min, z_max = min(pmol.df.z), max(pmol.df.z)
    z_center = z_min + abs(z_min-z_max)/2
    z_range = (abs(z_min-z_max)*size)
    return x_center, y_center, z_center, x_range, y_range, z_range

def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z):
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}')


dock_file(docker_command, example_receptor, example_active, *get_coordinates(example_crystal_ligand, 1))
