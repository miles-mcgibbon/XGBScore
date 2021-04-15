'''
This script loads binana.py output text files and extracts
the features to a .csv file
'''

import pandas as pd
import os
import io
import sys
import concurrent.futures
from tqdm import tqdm

# set global for multithreading
destination_path = None

def read_binana(binana_file): # load the text file and extract markdown table strings into a dictionary

    # read in binana output text file
    binana_file = open(binana_file, 'r')

    # cut out license and prose text
    binana_stats = binana_file.read().split('[- END INTRO -]')[1]

    # define headers for parsing
    file_headers = ['Atom-type pair counts within 2.5 angstroms:',
               'Atom-type pair counts within 4.0 angstroms:',
               'Ligand atom types:',
               'Summed electrostatic energy by atom-type pair, in J/mol:',
               'Number of rotatable bonds in the ligand:',
               'Active-site flexibility:',
               'Hydrogen bonds:',
               'Hydrophobic contacts (C-C):',
               'pi-pi stacking interactions:',
               'T-stacking (face-to-edge) interactions:',
               'Cation-pi interactions:',
               'Salt Bridges:']

    # define dictionary for populating
    tables = dict()

    # extract markdown tables from text file and add to tables dict
    for index, header in enumerate(file_headers):
        column_names = dict()
        if header == 'Number of rotatable bonds in the ligand:':
            table = binana_stats.split(header)[1].split(file_headers[index + 1])[0].strip().lstrip()
        elif index < (len(file_headers) - 1):
            table = binana_stats.split(header)[1].split(file_headers[index + 1])[0].strip().lstrip()
            table = table.split('Raw data')[0]
            table = io.StringIO(table)
            table = pd.read_table(table, sep="|")
            for column in list(table):
                column_names[column] = column.lstrip().rstrip()
                table = table.loc[~table[column].str.contains('----')]
            table.rename(columns = column_names, inplace = True)
        else:
            table = binana_stats.split(header)[1].strip().lstrip()
            table = table.split('Raw data')[0]
            table = io.StringIO(table)
            table = pd.read_table(table, sep="|")
            for column in list(table):
                column_names[column] = column.lstrip().rstrip()
                table = table.loc[~table[column].str.contains('----')]
            table.rename(columns = column_names, inplace = True)
        tables[header] = table

    return tables

def extract_features(file):

    # get markdown tables dictionary from file
    tables = read_binana(file)

    # define variables
    name = file.split('/')
    name = name[len(name)-1].replace('.txt','')

    # set headers to preserve column order in future uses
    df_headers = ['Name',
                '2.5 (HD, OA)', '2.5 (HD, HD)', '2.5 (C, HD)', '2.5 (HD, NA)', '2.5 (FE, HD)', '2.5 (A, HD)', '2.5 (HD, N)', '2.5 (OA, ZN)', '2.5 (N, ZN)', '2.5 (HD, ZN)', '2.5 (OA, OA)', '2.5 (CO, OA)', '2.5 (CO, HD)', '2.5 (C, C)', '2.5 (N, OA)', '2.5 (F, HD)', '2.5 (C, OA)', '2.5 (C, CL)', '2.5 (MG, OA)', '2.5 (NA, OA)', '2.5 (CA, OA)', '2.5 (CL, HD)', '2.5 (CL, MG)', '2.5 (A, OA)', '2.5 (C, MG)', '2.5 (C, N)', '2.5 (HD, SA)', '2.5 (MN, OA)', '2.5 (FE, NA)', '2.5 (FE, OA)', '2.5 (HD, P)', '2.5 (CA, HD)', '2.5 (HD, MN)', '2.5 (A, P)', '2.5 (F, OA)', '2.5 (C, ZN)', '2.5 (HD, MG)', '2.5 (N, N)', '2.5 (C, NA)', '2.5 (CO, NA)', '2.5 (A, C)', '2.5 (NA, ZN)', '2.5 (A, CL)', '2.5 (MG, N)', '2.5 (CA, NA)', '2.5 (MN, N)', '2.5 (CD, OA)', '2.5 (A, N)', '2.5 (CL, SA)', '2.5 (OA, SA)', '2.5 (HD, S)', '2.5 (BR, HD)', '2.5 (SA, SA)', '2.5 (SA, ZN)', '2.5 (N, SA)', '2.5 (NI, OA)', '2.5 (C, CA)', '2.5 (A, SA)', '2.5 (C, SA)', '2.5 (N, NA)', '2.5 (CL, N)', '2.5 (C, P)', '2.5 (A, A)', '2.5 (MN, NA)', '2.5 (CU, OA)', '2.5 (HD, NI)', '2.5 (C, CO)', '2.5 (A, MG)', '2.5 (CU, HD)', '2.5 (CO, N)', '2.5 (A, NA)', '2.5 (C, FE)', '2.5 (CL, OA)', '2.5 (N, NI)', '2.5 (BR, OA)', '2.5 (BR, C)', '2.5 (A, ZN)', '2.5 (C, F)', '2.5 (MG, P)', '2.5 (OA, P)', '2.5 (A, F)', '2.5 (MG, NA)', '2.5 (C, NI)', '2.5 (HD, K)', '2.5 (BR, FE)', '2.5 (A, BR)', '2.5 (C, MN)', '2.5 (F, N)', '2.5 (CU, NA)', '2.5 (CD, HD)', '2.5 (CU, N)', '2.5 (A, MN)', '2.5 (A, FE)', '2.5 (NA, NA)', '2.5 (F, MG)', '2.5 (A, CU)', '2.5 (CU, F)', '2.5 (FE, N)', '2.5 (C, CU)', '2.5 (NA, NI)', '2.5 (C, CD)', '2.5 (F, MN)', '2.5 (K, OA)', '2.5 (C, K)', '2.5 (S, ZN)', '2.5 (K, NA)', '2.5 (F, ZN)', '2.5 (CD, NA)', '2.5 (N, P)', '2.5 (F, SA)', '2.5 (F, FE)', '2.5 (OA, S)', '2.5 (BR, SA)', '2.5 (MG, SA)', '2.5 (BR, ZN)', '2.5 (A, NI)', '2.5 (CL, ZN)', '2.5 (CA, N)', '2.5 (MN, P)', '2.5 (HD, I)', '2.5 (P, ZN)', '2.5 (NA, SA)',
                '4.0 (C, C)', '4.0 (C, OA)', '4.0 (C, HD)', '4.0 (C, N)', '4.0 (N, OA)', '4.0 (HD, OA)', '4.0 (N, N)', '4.0 (HD, N)', '4.0 (HD, HD)', '4.0 (A, OA)', '4.0 (A, C)', '4.0 (A, A)', '4.0 (A, N)', '4.0 (A, HD)', '4.0 (C, SA)', '4.0 (HD, SA)', '4.0 (OA, OA)', '4.0 (C, ZN)', '4.0 (OA, ZN)', '4.0 (C, NA)', '4.0 (A, NA)', '4.0 (NA, OA)', '4.0 (N, NA)', '4.0 (HD, NA)', '4.0 (A, SA)', '4.0 (OA, SA)', '4.0 (BR, C)', '4.0 (BR, OA)', '4.0 (F, OA)', '4.0 (F, HD)', '4.0 (A, F)', '4.0 (F, N)', '4.0 (C, F)', '4.0 (A, FE)', '4.0 (C, FE)', '4.0 (FE, N)', '4.0 (FE, NA)', '4.0 (FE, OA)', '4.0 (C, P)', '4.0 (CL, N)', '4.0 (CL, HD)', '4.0 (C, CL)', '4.0 (CL, OA)', '4.0 (N, ZN)', '4.0 (HD, ZN)', '4.0 (N, SA)', '4.0 (SA, ZN)', '4.0 (OA, P)', '4.0 (N, P)', '4.0 (HD, P)', '4.0 (MG, OA)', '4.0 (A, ZN)', '4.0 (HD, S)', '4.0 (A, S)', '4.0 (N, S)', '4.0 (OA, S)', '4.0 (S, ZN)', '4.0 (A, P)', '4.0 (A, CL)', '4.0 (P, ZN)', '4.0 (C, S)', '4.0 (CO, P)', '4.0 (CO, HD)', '4.0 (F, ZN)', '4.0 (F, SA)', '4.0 (S, SA)', '4.0 (FE, HD)', '4.0 (NA, ZN)', '4.0 (CO, N)', '4.0 (A, CO)', '4.0 (C, CO)', '4.0 (CO, OA)', '4.0 (HD, MG)', '4.0 (C, MG)', '4.0 (CL, SA)', '4.0 (NA, SA)', '4.0 (SA, SA)', '4.0 (A, CA)', '4.0 (C, CA)', '4.0 (CA, N)', '4.0 (CA, OA)', '4.0 (CA, P)', '4.0 (CA, HD)', '4.0 (A, MG)', '4.0 (C, MN)', '4.0 (MN, N)', '4.0 (HD, MN)', '4.0 (MN, NA)', '4.0 (A, MN)', '4.0 (F, MN)', '4.0 (BR, N)', '4.0 (BR, HD)', '4.0 (A, BR)', '4.0 (MN, P)', '4.0 (MN, OA)', '4.0 (C, I)', '4.0 (I, SA)', '4.0 (A, I)', '4.0 (MG, NA)', '4.0 (MG, N)', '4.0 (MG, P)', '4.0 (HD, NI)', '4.0 (A, NI)', '4.0 (NI, OA)', '4.0 (I, OA)', '4.0 (C, NI)', '4.0 (MG, SA)', '4.0 (N, NI)', '4.0 (C, CD)', '4.0 (CD, N)', '4.0 (F, MG)', '4.0 (F, S)', '4.0 (C, CU)', '4.0 (P, SA)', '4.0 (BR, SA)', '4.0 (CL, ZN)', '4.0 (MN, SA)', '4.0 (BR, MG)', '4.0 (CA, F)', '4.0 (FE, P)', '4.0 (CL, CU)', '4.0 (CA, NA)', '4.0 (NA, NI)', '4.0 (CU, N)', '4.0 (CU, HD)', '4.0 (CU, OA)', '4.0 (MG, S)', '4.0 (A, CU)', '4.0 (NA, P)', '4.0 (CO, NA)', '4.0 (CU, NA)', '4.0 (BR, ZN)', '4.0 (HD, I)', '4.0 (I, N)', '4.0 (F, NA)', '4.0 (NA, NA)', '4.0 (FE, SA)', '4.0 (C, K)', '4.0 (K, OA)', '4.0 (K, N)', '4.0 (HD, K)', '4.0 (K, P)', '4.0 (FE, S)', '4.0 (CD, OA)', '4.0 (CD, HD)', '4.0 (CA, SA)', '4.0 (A, HG)', '4.0 (HG, SA)', '4.0 (CD, P)', '4.0 (CO, S)', '4.0 (CA, S)', '4.0 (CL, MG)', '4.0 (CU, F)', '4.0 (NI, SA)', '4.0 (MN, S)', '4.0 (BR, NA)', '4.0 (F, FE)', '4.0 (F, K)', '4.0 (K, SA)', '4.0 (CL, FE)', '4.0 (A, CD)', '4.0 (CL, MN)', '4.0 (NA, S)', '4.0 (CA, CL)', '4.0 (BR, FE)', '4.0 (A, K)', '4.0 (K, NA)', '4.0 (CD, NA)', '4.0 (CL, NA)', '4.0 (CO, SA)', '4.0 (HG, OA)', '4.0 (CL, NI)', '4.0 (CD, SA)', '4.0 (CD, S)', '4.0 (BR, MN)',
                'LA C', 'LA N', 'LA HD', 'LA OA', 'LA A', 'LA SA', 'LA NA', 'LA BR', 'LA F', 'LA P', 'LA CL', 'LA S', 'LA I',
                'ElSum (C, C)', 'ElSum (C, OA)', 'ElSum (C, HD)', 'ElSum (C, N)', 'ElSum (N, OA)', 'ElSum (HD, OA)', 'ElSum (N, N)', 'ElSum (HD, N)', 'ElSum (HD, HD)', 'ElSum (A, OA)', 'ElSum (A, C)', 'ElSum (A, A)', 'ElSum (A, N)', 'ElSum (A, HD)', 'ElSum (C, SA)', 'ElSum (HD, SA)', 'ElSum (OA, OA)', 'ElSum (C, ZN)', 'ElSum (OA, ZN)', 'ElSum (C, NA)', 'ElSum (A, NA)', 'ElSum (NA, OA)', 'ElSum (N, NA)', 'ElSum (HD, NA)', 'ElSum (A, SA)', 'ElSum (OA, SA)', 'ElSum (BR, C)', 'ElSum (BR, OA)', 'ElSum (F, OA)', 'ElSum (F, HD)', 'ElSum (A, F)', 'ElSum (F, N)', 'ElSum (C, F)', 'ElSum (A, FE)', 'ElSum (C, FE)', 'ElSum (FE, N)', 'ElSum (FE, HD)', 'ElSum (FE, NA)', 'ElSum (FE, OA)', 'ElSum (C, P)', 'ElSum (CL, N)', 'ElSum (CL, HD)', 'ElSum (C, CL)', 'ElSum (CL, OA)', 'ElSum (N, ZN)', 'ElSum (HD, ZN)', 'ElSum (N, SA)', 'ElSum (SA, ZN)', 'ElSum (OA, P)', 'ElSum (N, P)', 'ElSum (HD, P)', 'ElSum (MG, OA)', 'ElSum (A, ZN)', 'ElSum (HD, S)', 'ElSum (A, S)', 'ElSum (N, S)', 'ElSum (OA, S)', 'ElSum (S, ZN)', 'ElSum (A, P)', 'ElSum (A, CL)', 'ElSum (P, ZN)', 'ElSum (C, S)', 'ElSum (CO, P)', 'ElSum (CO, OA)', 'ElSum (CO, HD)', 'ElSum (F, ZN)', 'ElSum (F, SA)', 'ElSum (S, SA)', 'ElSum (NA, ZN)', 'ElSum (CO, N)', 'ElSum (A, CO)', 'ElSum (C, CO)', 'ElSum (HD, MG)', 'ElSum (C, MG)', 'ElSum (CL, SA)', 'ElSum (NA, SA)', 'ElSum (SA, SA)', 'ElSum (A, CA)', 'ElSum (C, CA)', 'ElSum (CA, N)', 'ElSum (CA, OA)', 'ElSum (CA, P)', 'ElSum (CA, HD)', 'ElSum (A, MG)', 'ElSum (CL, MG)', 'ElSum (C, MN)', 'ElSum (MN, N)', 'ElSum (MN, OA)', 'ElSum (HD, MN)', 'ElSum (MN, NA)', 'ElSum (A, MN)', 'ElSum (F, MN)', 'ElSum (BR, N)', 'ElSum (BR, HD)', 'ElSum (A, BR)', 'ElSum (MN, P)', 'ElSum (C, I)', 'ElSum (I, SA)', 'ElSum (A, I)', 'ElSum (MG, NA)', 'ElSum (MG, N)', 'ElSum (MG, P)', 'ElSum (HD, NI)', 'ElSum (CO, NA)', 'ElSum (A, NI)', 'ElSum (NI, OA)', 'ElSum (I, OA)', 'ElSum (C, NI)', 'ElSum (MG, SA)', 'ElSum (CA, NA)', 'ElSum (N, NI)', 'ElSum (C, CD)', 'ElSum (CD, N)', 'ElSum (CD, OA)', 'ElSum (F, MG)', 'ElSum (F, S)', 'ElSum (C, CU)', 'ElSum (P, SA)', 'ElSum (BR, SA)', 'ElSum (CL, ZN)', 'ElSum (MN, SA)', 'ElSum (BR, MG)', 'ElSum (CA, F)', 'ElSum (FE, P)', 'ElSum (CL, CU)', 'ElSum (NA, NI)', 'ElSum (CU, N)', 'ElSum (CU, OA)', 'ElSum (CU, HD)', 'ElSum (MG, S)', 'ElSum (A, CU)', 'ElSum (NA, P)', 'ElSum (CU, NA)', 'ElSum (BR, ZN)', 'ElSum (HD, I)', 'ElSum (I, N)', 'ElSum (F, NA)', 'ElSum (NA, NA)', 'ElSum (FE, SA)', 'ElSum (C, K)', 'ElSum (K, OA)', 'ElSum (K, N)', 'ElSum (HD, K)', 'ElSum (K, P)', 'ElSum (FE, S)', 'ElSum (CD, HD)', 'ElSum (CA, SA)', 'ElSum (BR, FE)', 'ElSum (A, HG)', 'ElSum (HG, SA)', 'ElSum (CD, P)', 'ElSum (CO, S)', 'ElSum (CA, S)', 'ElSum (CU, F)', 'ElSum (NI, SA)', 'ElSum (MN, S)', 'ElSum (BR, NA)', 'ElSum (F, FE)', 'ElSum (F, K)', 'ElSum (K, SA)', 'ElSum (CL, FE)', 'ElSum (A, CD)', 'ElSum (CL, MN)', 'ElSum (NA, S)', 'ElSum (CA, CL)', 'ElSum (A, K)', 'ElSum (K, NA)', 'ElSum (CD, NA)', 'ElSum (CL, NA)', 'ElSum (CO, SA)', 'ElSum (HG, OA)', 'ElSum (CL, NI)', 'ElSum (CD, SA)', 'ElSum (CD, S)', 'ElSum (BR, MN)',
                'BPF ALPHA SIDECHAIN', 'BPF ALPHA BACKBONE', 'BPF BETA SIDECHAIN', 'BPF BETA BACKBONE', 'BPF OTHER SIDECHAIN', 'BPF OTHER BACKBONE',
                'HC ALPHA SIDECHAIN', 'HC ALPHA BACKBONE', 'HC BETA SIDECHAIN', 'HC BETA BACKBONE', 'HC OTHER SIDECHAIN', 'HC OTHER BACKBONE',
                'HB ALPHA SIDECHAIN LIGAND', 'HB ALPHA BACKBONE LIGAND', 'HB BETA SIDECHAIN LIGAND', 'HB BETA BACKBONE LIGAND', 'HB OTHER SIDECHAIN LIGAND', 'HB OTHER BACKBONE LIGAND',
                'HB ALPHA SIDECHAIN RECEPTOR', 'HB ALPHA BACKBONE RECEPTOR', 'HB BETA SIDECHAIN RECEPTOR', 'HB BETA BACKBONE RECEPTOR', 'HB OTHER SIDECHAIN RECEPTOR', 'HB OTHER BACKBONE RECEPTOR',
                'SB ALPHA', 'SB BETA', 'SB OTHER',
                'piStack ALPHA', 'piStack BETA', 'piStack OTHER',
                'tStack ALPHA', 'tStack BETA', 'tStack OTHER',
                'catPi ALPHA LIGAND', 'catPi BETA LIGAND', 'catPi OTHER LIGAND',
                'catPi ALPHA RECEPTOR', 'catPi BETA RECEPTOR', 'catPi OTHER RECEPTOR',
                'nRot']

    # make empty dataframe for populating with forced headers
    binana_data = pd.DataFrame(columns = df_headers)

    # set up dictionary to populate with binana features for structure
    data = {'Name': str(name)}

    # iterate through markdown tables extracting each category of features
    for key, value in tables.items():

        if key == 'Atom-type pair counts within 2.5 angstroms:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'2.5 ({atom1}, {atom2})'] = tally


        elif key == 'Atom-type pair counts within 4.0 angstroms:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'4.0 ({atom1}, {atom2})'] = tally

        elif key == 'Ligand atom types:':
            for index in value.index:
                atom = (value.iloc[index-1,0]).lstrip().rstrip()
                data[f'LA {atom}'] = 1

        elif key == 'Summed electrostatic energy by atom-type pair, in J/mol:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = float(value.iloc[index-1,2])
                data[f'ElSum ({atom1}, {atom2})'] = tally

        elif key == 'Active-site flexibility:':
            for index in value.index:
                chain = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'BPF {struc} {chain}'] = tally

        elif key == 'Hydrogen bonds:':
            for index in value.index:
                donor = (value.iloc[index-1,0]).lstrip().rstrip()
                chain = (value.iloc[index-1,1]).lstrip().rstrip()
                struc = (value.iloc[index-1,2]).lstrip().rstrip()
                tally = int(value.iloc[index-1,3])
                data[f'HB {struc} {chain} {donor}'] = tally

        elif key == 'Hydrophobic contacts (C-C):':
            for index in value.index:
                chain = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'HC {struc} {chain}'] = tally

        elif key == 'pi-pi stacking interactions:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'piStack {struc}'] = tally

        elif key == 'T-stacking (face-to-edge) interactions:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'tStack {struc}'] = tally

        elif key == 'Cation-pi interactions:':
            for index in value.index:
                char = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'catPi {struc} {char}'] = tally

        elif key == 'Salt Bridges:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'SB {struc}'] = tally

        elif key == 'Number of rotatable bonds in the ligand:':
            data['nRot'] = value

    # add data to dataframe
    binana_data = binana_data.append(data, ignore_index = True)

    # fill absent features with zeroes
    binana_data = binana_data.fillna(0)

    # print any unexpected columns (this should not happen)
    if len(binana_data.columns) != 526:
        print('Col Num')
        print(len(binana_data.columns))
        print(name)

    # save the dataframe with added row
    if os.path.isfile(destination_path):
        binana_data.to_csv(destination_path, mode = 'a', index = False, header = False)
    else:
        binana_data.to_csv(destination_path, index = False)

def parse_args(args): # parse CLI user inputs

    binana_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-out') + 1]

    return binana_files_path, destination_path

def main(): # run script using CLI

    global destination_path

    binana_files_path, destination_path = parse_args(sys.argv)

    # get list of absolute filepaths
    binana_files = [(binana_files_path + text_file) for text_file in os.listdir(binana_files_path)]

    # extract binana features for each file and save
    for file in tqdm(binana_files):
        binana_data = extract_features(file)

if __name__ == '__main__':
    main()
