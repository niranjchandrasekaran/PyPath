from parameter.constants import Constant
from scipy.spatial.distance import pdist, squareform
import numpy as np
from parameter.mass import Mass
from parameter.sidechain import SideChain


class BuildHessian(object):

    def __init__(self):
        pass

    def hessian(self, coord, pdb, calpha):
        if calpha:
            hessian = self.hessian_amweh(coord, pdb.natoms, pdb.aaname)
        else:
            hessian = self.hessian_all_atom(coord, pdb)

        return hessian

    def hessian_all_atom(self, coord, pdb):
        anm = self.hessian_anm(coord, pdb.natoms)
        torsion_backbone = self.hessian_backbone_torsional(coord, pdb.atname, pdb.natoms)
        torsion_sidechain = self.hessian_sidechain_torsional(coord, pdb.atname, pdb.aaname, pdb.aano, pdb.natoms)

        hessian = []

        for i in range(len(anm)):
            hessian.append([])
            for j in range(len(anm)):
                hessian[i].append(anm[i][j] + torsion_backbone[i][j] + torsion_sidechain[i][j])

        return hessian

    def hessian_anm(self, coord, natoms):
        constant = Constant()
        hessian = [[0 for y in range(constant.dim * natoms)] for x in range(constant.dim * natoms)]
        distance = squareform(pdist(np.asarray(coord), 'euclidean'))

        for atom_i in range(natoms):
            for atom_j in range(natoms):
                if atom_i == atom_j:
                    continue
                else:
                    if distance[atom_i][atom_j] <= constant.cutoff:
                        for x in range(constant.dim):
                            for y in range(constant.dim):
                                hessian[(3 * atom_i) + x][(3 * atom_j) + y] = -constant.force_constant * (
                                        coord[atom_i][x] - coord[atom_j][x]) * (coord[atom_i][y] - coord[atom_j][
                                    y]) / distance[atom_i][atom_j]
                                hessian[(3 * atom_i) + x][(3 * atom_i) + y] += -hessian[(3 * atom_i) + x][
                                    (3 * atom_j) + y]

        return hessian

    def hessian_backbone_torsional(self, coord, atom_name, natoms):
        constant = Constant()
        backbone = ['N', 'CA', 'C']
        hessian = [[0 for y in range(constant.dim * natoms)] for x in range(constant.dim * natoms)]

        backbone = [atom in backbone for atom in atom_name]

        backbone_torsion_atoms = []

        for bb_i in range(len(backbone)):
            if backbone[bb_i] and atom_name[bb_i] != 'CA':
                count = 0
                backbone_torsion_atoms.append([])
                for bb_j in range(bb_i, len(backbone)):
                    if backbone[bb_j]:
                        backbone_torsion_atoms[-1].append(bb_j)
                        count += 1
                    if count == 4:
                        break
                if len(backbone_torsion_atoms[-1]) < 4:
                    backbone_torsion_atoms.remove(backbone_torsion_atoms[-1])

        for tetrad in range(len(backbone_torsion_atoms)):
            tetrad_coord = []
            for atom in range(len(backbone_torsion_atoms[tetrad])):
                tetrad_coord.append(coord[backbone_torsion_atoms[tetrad][atom]])

            hess = self.torsional_hessian(tetrad_coord)

            for atom_i in range(0, 4):
                for atom_j in range(0, 4):
                    for i in range(0, 3):
                        for j in range(0, 3):
                            hessian[(backbone_torsion_atoms[tetrad][atom_i] * 3) + i][
                                (backbone_torsion_atoms[tetrad][atom_j] * 3) + j] += hess[3 * atom_i + i][
                                3 * atom_j + j]

        return hessian

    def hessian_sidechain_torsional(self, coord, atom_name, amino_acid, amino_acid_no, natoms):
        constant = Constant()
        sidechain = SideChain()
        amino_acid_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO',
                           'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        hessian = [[0 for y in range(constant.dim * natoms)] for x in range(constant.dim * natoms)]

        amino_acid_present = [aa in amino_acid_list for aa in amino_acid]

        start = 0

        sidechain_torsion_atoms = []

        for aa in range(len(amino_acid_no)-1):
            if amino_acid_present[aa]:
                if aa == len(amino_acid_no)-2:
                    end = len(amino_acid_no)
                    if amino_acid[aa] != 'GLY' and amino_acid[aa] != 'ALA':
                        tetrad_atoms_list = sidechain.amino_acid(amino_acid[aa], atom_name, start, end)
                        for tetrad in range(len(tetrad_atoms_list)):
                            sidechain_torsion_atoms.append(tetrad_atoms_list[tetrad])
                else:
                    if amino_acid_no[aa+1] != amino_acid_no[aa]:
                        end = aa+1
                        if amino_acid[aa] != 'GLY' and amino_acid[aa] != 'ALA':
                            tetrad_atoms_list = sidechain.amino_acid(amino_acid[aa], atom_name, start, end)
                            for tetrad in range(len(tetrad_atoms_list)):
                                sidechain_torsion_atoms.append(tetrad_atoms_list[tetrad])
                        start = aa+1

        for tetrad in range(len(sidechain_torsion_atoms)):
            tetrad_coord = []
            for atom in range(len(sidechain_torsion_atoms[tetrad])):
                tetrad_coord.append(coord[sidechain_torsion_atoms[tetrad][atom]])

            hess = self.torsional_hessian(tetrad_coord)

            for atom_i in range(0, 4):
                for atom_j in range(0, 4):
                    for i in range(0, 3):
                        for j in range(0, 3):
                            hessian[(sidechain_torsion_atoms[tetrad][atom_i] * 3) + i][
                                (sidechain_torsion_atoms[tetrad][atom_j] * 3) + j] += hess[3 * atom_i + i][
                                3 * atom_j + j]

        return hessian

    def torsional_hessian(self, coord):
        x1, y1, z1 = coord[0]
        x2, y2, z2 = coord[1]
        x3, y3, z3 = coord[2]
        x4, y4, z4 = coord[3]

        a, b, c = [[] for x in range(3)]

        for i in range(3):
            a.append(coord[1][i] - coord[0][i])
            b.append(coord[2][i] - coord[1][i])
            c.append(coord[3][i] - coord[2][1])

        v1 = np.cross(a, b)
        v2 = np.cross(b, c)

        v1m = np.linalg.norm(v1)
        v2m = np.linalg.norm(v2)

        r = float(np.dot(v1, v2)) / float(v1m * v2m)

        R = 0.00030461741 * float(1) / float(((v1m * v2m) ** 2) * np.sqrt(1 - r ** 2))

        d = []

        k1 = (y2 - y1) * (z3 - z2) - (y3 - y2) * (z2 - z1)
        k2 = (x2 - x1) * (z3 * z2) - (x3 - x2) * (z2 - z1)
        k3 = (x2 - x1) * (y3 - y2) - (x3 - x2) * (y2 - y1)

        l1 = (y3 - y2) * (z4 - z3) - (z3 - z2) * (y4 - y3)
        l2 = (x3 - x2) * (z4 - z3) - (z3 - z2) * (x4 - x3)
        l3 = (x3 - x2) * (y4 - y3) - (y3 - y2) * (x4 - x3)

        d.append(R * (v1m * v2m * (l2 * (z2 - z3) + l3 * (y2 - y3)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k2 * (z2 - z3) + k3 * (y2 - y3)))))
        d.append(R * (v1m * v2m * (l1 * (z2 - z3) + l3 * (x3 - x2)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (z2 - z3) + k3 * (x3 - x2)))))
        d.append(R * (v1m * v2m * (l1 * (y3 - y2) + l2 * (x3 - x2)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (y3 - y2) + k2 * (x3 - x2)))))

        d.append(R * (
                v1m * v2m * (l2 * (z3 - z1) + l3 * (y3 - y1) + k2 * (z3 - z4) + k3 * (y3 - y4)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (z3 - z1) + k2 * (y3 - y1)) + float(v1m) / float(v2m) * (
                l2 * (z3 - z4) + l3 * (y3 - y4)))))
        d.append(R * (
                v1m * v2m * (l1 * (z3 - z1) + l3 * (x1 - x3) + k1 * (z3 - z4) + k3 * (x4 - x3)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (z3 - z1) + k3 * (x1 - x3)) + float(v1m) / float(v2m) * (
                l1 * (z3 - z4) + l3 * (x4 - x3)))))
        d.append(R * (
                v1m * v2m * (l1 * (y1 - y3) + l2 * (x1 - x3) + k1 * (y4 - y3) + k2 * (x4 - x3)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (y1 - y3) + k2 * (x1 - x3)) + float(v1m) / float(v2m) * (
                l1 * (y4 - y3) + l2 * (x4 - x3)))))

        d.append(R * (
                v1m * v2m * (l2 * (z1 - z2) + l3 * (y1 - y2) + k2 * (z4 - z2) + k3 * (y4 - y2)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k2 * (z1 - z2) + k3 * (y1 - y2)) + float(v1m) / float(v2m) * (
                l2 * (z4 - z2) + l3 * (y4 - y2)))))
        d.append(R * (
                v1m * v2m * (l1 * (z1 - z2) + l3 * (x2 - x1) + k1 * (z4 - z2) + k3 * (x2 - x4)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (z1 - z2) + k3 * (x2 - x1)) + float(v1m) / float(v2m) * (
                l1 * (z4 - z2) + l3 * (x2 - x4)))))
        d.append(R * (
                v1m * v2m * (l1 * (y2 - y1) + l2 * (x2 - x1) + k1 * (y2 - y4) + k2 * (x2 - x4)) - np.dot(v1, v2) * (
                float(v2m) / float(v1m) * (k1 * (y2 - y1) + k2 * (x2 - x1)) + float(v1m) / float(v2m) * (
                l1 * (y2 - y4) + l2 * (x2 - x4)))))

        d.append(R * (v1m * v2m * (k2 * (z2 - z3) + k3 * (y2 - y3)) - np.dot(v1, v2) * (
                float(v1m) / float(v2m) * (l2 * (z2 - z3) + l3 * (y2 - y3)))))
        d.append(R * (v1m * v2m * (k1 * (z2 - z3) + k3 * (x3 - x2)) - np.dot(v1, v2) * (
                float(v1m) / float(v2m) * (l1 * (z2 - z3) + l3 * (x3 - x2)))))
        d.append(R * (v1m * v2m * (k1 * (y3 - y2) + k2 * (x3 - x2)) - np.dot(v1, v2) * (
                float(v1m) / float(v2m) * (l1 * (y3 - y2) + l2 * (x3 - x2)))))

        hess = [[d[i] * d[j] for j in range(len(d))] for i in range(len(d))]

        return hess

    def hessian_amweh(self, coord, natoms, amino_acid):
        constant = Constant()
        mass = Mass()
        hessian = [[0 for y in range(constant.dim * natoms)] for x in range(constant.dim * natoms)]
        distance = squareform(pdist(np.asarray(coord), 'euclidean'))

        for atom_i in range(natoms):
            for atom_j in range(natoms):
                if atom_i == atom_j:
                    continue
                else:
                    if distance[atom_i][atom_j] <= 4:
                        force_constant = 8.6e2 * distance[atom_i][atom_j] - 2.39e3
                    else:
                        force_constant = 128e4 / distance[atom_i][atom_j] ** 6
                    for x in range(constant.dim):
                        for y in range(constant.dim):
                            hessian[(3 * atom_i) + x][(3 * atom_j) + y] = -(
                                    force_constant * (coord[atom_i][x] - coord[atom_j][y]) * (
                                    coord[atom_i][y] - coord[atom_j][y])) / (distance[atom_i][atom_j] * (
                                    np.sqrt(mass.amino_acid_mass(amino_acid[atom_i])) * np.sqrt(
                                mass.amino_acid_mass(amino_acid[atom_j]))))

        return hessian
