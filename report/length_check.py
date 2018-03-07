import sys

class LengthCheck(object):

    def __init__(self, len_pdb1, len_pdb2):
        if len_pdb1 != len_pdb2:
            print('\n@> The two molecules are unequal in size.\n')
            sys.exit(1)
        else:
            if len_pdb1 > 1:
                print('\n@> There are %d atom(s) in your molecule.\n' % len_pdb1)
            else:
                print('\n@> Path requires the molecules to be at least diatomic.\n' )
                sys.exit(1)