class PDBWrite(object):

    def __init__(self, coord, pdb, filename):
        fout = open(filename, 'w')

        for atom in range(pdb.natoms):
            fout.write('ATOM   %4d %-4s %s %s%4d    %8.3f%8.3f%8.3f\n' % (
            (atom + 1), pdb.atname[atom], pdb.aaname[atom], pdb.chain[atom], int(pdb.aano[atom]), float(coord[atom][0]),
            float(coord[atom][1]), float(coord[atom][2])))

        fout.write('TER\n')
        fout.write('ENDMDL\n')
        fout.close()


class PDBTrajectoryWrite(object):

    def __init__(self):
        pass
