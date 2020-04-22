import sys

class LengthCheck(object):

    def __init__(self, len_pdb1, len_pdb2):
        """
        Checks to make sure the number of atoms of the two end states are equal

        :param len_pdb1: number of atoms of the initial state
        :param len_pdb2: number of atoms of the final state
        """
        if len_pdb1 != len_pdb2:
            print('@> The two molecules are unequal in size.\n')
            sys.exit(1)
        else:
            if len_pdb1 > 1:
                print('@> There are %d atom(s) in your molecule.\n' % len_pdb1)
            else:
                print('@> Path requires the molecules to be at least diatomic.\n' )
                sys.exit(1)