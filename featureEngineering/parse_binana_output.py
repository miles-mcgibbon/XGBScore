import pandas as pd
import io

binana_file = open('/home/milesm/Desktop/7cpa/7cpa_BINANA_Data.txt', 'r')

binana_stats = binana_file.read().split('[- END INTRO -]')[1]

headers = ['Atom-type pair counts within 2.5 angstroms:',
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

tables = {}
for index, header in enumerate(headers):
    if header == 'Number of rotatable bonds in the ligand:':
        table = binana_stats.split(header)[1].split(headers[index + 1])[0].strip().lstrip()
    elif index < (len(headers) - 1):
        table = binana_stats.split(header)[1].split(headers[index + 1])[0].strip().lstrip()
        table = table.split('Raw data')[0]
        table = io.StringIO(table)
        table = pd.read_table(table, sep="|")
        for column in list(table):
            table = table.loc[~table[column].str.contains('----')]
    else:
        table = binana_stats.split(header)[1].strip().lstrip()
        table = table.split('Raw data')[0]
        table = io.StringIO(table)
        table = pd.read_table(table, sep="|")
        for column in list(table):
            table = table.loc[~table[column].str.contains('----')]
    tables[header] = table
