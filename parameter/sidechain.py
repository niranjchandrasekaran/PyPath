class SideChain(object):

    def __init__(self):
        pass

    def amino_acid(self, amino_acid, atom_list, start, end):
        if amino_acid == 'CYS':
            tetrads = [('N', 'CA', 'CB', 'SG')]
        elif amino_acid == 'ASP':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')]
        elif amino_acid == 'GLU':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')]
        elif amino_acid == 'PHE':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')]
        elif amino_acid == 'HIS':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'ND1')]
        elif amino_acid == 'ILE':
            tetrads = [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')]
        elif amino_acid == 'LYS':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'CE'),
                      ('CG', 'CD', 'CE', 'NZ')]
        elif amino_acid == 'LEU':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')]
        elif amino_acid == 'MET':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'SD'), ('CB', 'CG', 'SD', 'CE')]
        elif amino_acid == 'ASN':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')]
        elif amino_acid == 'PRO':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')]
        elif amino_acid == 'GLN':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'OE1')]
        elif amino_acid == 'ARG':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD'), ('CB', 'CG', 'CD', 'NE'),
                      ('CG', 'CD', 'NE', 'CZ'), ('CD', 'NE', 'CZ', 'NH1')]
        elif amino_acid == 'SER':
            tetrads = [('N', 'CA', 'CB', 'OG')]
        elif amino_acid == 'THR':
            tetrads = [('N', 'CA', 'CB', 'OG1')]
        elif amino_acid == 'VAL':
            tetrads = [('N', 'CA', 'CB', 'CG1')]
        elif amino_acid == 'TRP':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')]
        elif amino_acid == 'TYR':
            tetrads = [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')]

        tetrad_atom_list = []

        for tetrad in range(len(tetrads)):
            tetrad_atom_list.append([])
            for atom in range(len(tetrads[tetrad])):
                for atom_id in range(start, end):
                    if tetrads[tetrad][atom] == atom_list[atom_id]:
                        tetrad_atom_list[tetrad].append(atom_id)

        return tetrad_atom_list
