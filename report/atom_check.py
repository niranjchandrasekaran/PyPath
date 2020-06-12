import sys
import pandas as pd


class AtomCheck(object):

    def __init__(self, pdb, calpha):
        """
        Checks to make sure that all the necessary atom are present in the PDB

        :param PDB file to check
        :calpha C-alpha flag
        """

        pdb_dict = {
            'chain': pdb.chain,
            'aano': pdb.aano,
            'aaname': pdb.aaname,
            'atname': pdb.atname
        }

        if not calpha:
            pdb_df = (
                pd.DataFrame(pdb_dict)
                    .groupby(['chain', 'aano', 'aaname']).atname.apply(list)
                    .reset_index()
                    .assign(expected=lambda x: x.aaname, presence=True)
            )

            pdb_df.expected = pdb_df.expected.apply(self.atomlist)
            pdb_df.presence = self.check_presence(list(pdb_df.atname), list(pdb_df.expected))

            missing_atoms = pdb_df.loc[~pdb_df.presence].reset_index(drop=True)

            if missing_atoms.shape[0] > 0:
                print(
                    '@> There are sidechain atoms missing in some amino acids. A list of these amino acids have been'
                    ' written to the file missing_sidechain_atoms.txt\n')

                self.write_missing(missing_atoms)
                sys.exit(1)
            else:
                print(f'@> All the required atoms are present in {pdb.name}.pdb\n')

    @staticmethod
    def atomlist(amino_acid):
        atom_list = []
        if amino_acid == 'CYS':
            atom_list = ['N', 'CA', 'CB', 'SG']
        elif amino_acid == 'ASP':
            atom_list = ['N', 'CA', 'CB', 'CG', 'OD1']
        elif amino_acid == 'GLU':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD', 'OE1']
        elif amino_acid == 'PHE':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD1']
        elif amino_acid == 'HIS':
            atom_list = ['N', 'CA', 'CB', 'CG', 'ND1']
        elif amino_acid == 'ILE':
            atom_list = ['N', 'CA', 'CB', 'CG1', 'CD1']
        elif amino_acid == 'LYS':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ']
        elif amino_acid == 'LEU':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD1']
        elif amino_acid == 'MET':
            atom_list = ['N', 'CA', 'CB', 'CG', 'SD', 'CE']
        elif amino_acid == 'ASN':
            atom_list = ['N', 'CA', 'CB', 'CG', 'OD1']
        elif amino_acid == 'PRO':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD']
        elif amino_acid == 'GLN':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD', 'OE1']
        elif amino_acid == 'ARG':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1']
        elif amino_acid == 'SER':
            atom_list = ['N', 'CA', 'CB', 'OG']
        elif amino_acid == 'THR':
            atom_list = ['N', 'CA', 'CB', 'OG1']
        elif amino_acid == 'VAL':
            atom_list = ['N', 'CA', 'CB', 'CG1']
        elif amino_acid == 'TRP':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD1']
        elif amino_acid == 'TYR':
            atom_list = ['N', 'CA', 'CB', 'CG', 'CD1']

        return atom_list

    @staticmethod
    def check_presence(atname, expected):
        presence_array = [set(expected[_]).issubset(set(atname[_])) for _ in range(len(atname))]
        return presence_array

    @staticmethod
    def write_missing(df):
        fout = open('missing_sidechain_atoms.txt', 'w')

        for row in range(len(df)):
            fout.write(f'{df.loc[row].aaname} {df.loc[row].aano} in chain {df.loc[row].chain} has the atoms '
                       f'{df.loc[row].atname} while the atoms {df.loc[row].expected} are expected.\n')

        fout.close()

