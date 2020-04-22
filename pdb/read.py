class PDBRead(object):

    def __init__(self, filename, calpha):
        """
        Generate the pdb object

        :param filename: name of the PDB file
        :param calpha: True for CA only simulation and False for all atom simulation
        """
        with open(filename, 'r') as fopen:
            data = [self.read_each_line(line) for line in fopen if line.startswith('ATOM') or line.startswith('HETATM')]
        if calpha:
            pdb = [data[line] for line in range(len(data)) if data[line][2] == 'CA']
        else:
            pdb = data
        self.name = filename[:filename.index('.pdb')].split('/')[-1]
        self.atno = [int(pdb[atom][1]) for atom in range(len(pdb))]
        self.atname = [(pdb[atom][2].rstrip()).lstrip() for atom in range(len(pdb))]
        self.aaname = [pdb[atom][3].rstrip() for atom in range(len(pdb))]
        self.chain = [pdb[atom][4] for atom in range(len(pdb))]
        self.aano = [int(pdb[atom][5]) for atom in range(len(pdb))]
        self.coord = [[float(pdb[atom][6]), float(pdb[atom][7]), float(pdb[atom][8])] for atom in range(len(pdb))]
        self.natoms = len(self.atno)

    def read_each_line(self, line):
        """
        :param line: Each line of a PDB file
        :return: split_line: Columns of PDB separated by commas
        """
        split_line = [line[0:6].rstrip(), line[6:11].lstrip(), line[12:16].rstrip().lstrip(), line[17:20].lstrip(),
                      line[21:22].lstrip(), line[22:26].lstrip(), line[30:38].lstrip(), line[38:46].lstrip(),
                      line[46:54].lstrip()]

        return split_line
