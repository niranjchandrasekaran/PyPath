class PDBRead(object):

    def __init__(self, filename, calpha):
        with open(filename, 'r') as fopen:
            data = [line.rstrip().split() for line in fopen if
                   not line.rstrip() == 'TER' and not line.rstrip() == 'ENDMDL']

        if calpha:
            pdb = [data[line] for line in range(len(data)) if (data[line][2].rstrip()).lstrip() == 'CA']
        else:
            pdb = data
        self.name = filename[:filename.index('.pdb')].split('/')[-1]
        self.atno = [int(pdb[atom][1]) for atom in range(len(pdb))]
        self.atname = [pdb[atom][2] for atom in range(len(pdb))]
        self.aaname = [pdb[atom][3] for atom in range(len(pdb))]
        self.chain = [pdb[atom][4] for atom in range(len(pdb))]
        self.aano = [int(pdb[atom][5]) for atom in range(len(pdb))]
        self.coord = [[float(pdb[atom][6]), float(pdb[atom][7]), float(pdb[atom][8])] for atom in range(len(pdb))]
        self.natoms = len(self.atno)