class Mass(object):

    def __init__(self):
        pass

    def amino_acid_mass(self, amino_acid):
        mass_list = {'ALA': 71.00788,
                     'CYS': 103.1388,
                     'ASP': 115.0886,
                     'GLU': 129.1155,
                     'PHE': 147.1766,
                     'GLY': 57.0519,
                     'HIS': 137.1411,
                     'ILE': 113.1594,
                     'LYS': 128.1741,
                     'LEU': 113.1594,
                     'MET': 131.1926,
                     'ASN': 114.1038,
                     'PRO': 97.1167,
                     'GLN': 128.1307,
                     'ARG': 156.1875,
                     'SER': 87.0782,
                     'THR': 101.1051,
                     'VAL': 99.1326,
                     'TRP': 186.2132,
                     'TYR': 163.1760,
                     'TRX': 186.2132}

        for mass in mass_list:
            if amino_acid == mass:
                return mass_list[mass]

    def atom_mass(self, atom):
        mass_list = {'C': 12.08, 'N': 14.01, 'O': 15.99, 'M': 24.31, 'P': 30.97, 'S': 32.06, 'Z': 65.38}

        for mass in mass_list:
            if mass == atom:
                return mass_list[mass]
