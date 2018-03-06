import numpy as np


class ThermoDynamics(object):

    def ___init__(self):
        pass

    def energy(self, hessian, work):
        return 0.5*np.dot(work.T, np.dot(np.asarray(hessian), work))

    def work(self, coord_1, coord_2):
        return np.asarray(coord_1).ravel() - np.asarray(coord_2).ravel()