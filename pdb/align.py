from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

class Align(object):

    def __init__(self, static, moving):
        sup = SVDSuperimposer()
        sup.set(np.asarray(static), np.asarray(moving))
        sup.run()

        rot, trans = sup.get_rotran()

        self.rms = sup.get_rms()

        self.static = static

        self.moving = [np.dot(np.asarray(moving[atom]), rot) + trans for atom in range(len(moving))]