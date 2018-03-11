import numpy as np


class ThermoDynamics(object):

    def ___init__(self):
        """
        Calculates the thermodynamics of the structure
        """
        pass

    def energy(self, hessian, work):
        """
        Calculate the energy relative to another

        :param hessian: hessian of the end state
        :param work: work vector of structure relative to the end state
        :return: returns the relative energy
        """
        return 0.5*np.dot(work.T, np.dot(np.asarray(hessian), work))

    def work(self, coord_1, coord_2):
        """
        Calculates the work vector relative to another

        :param coord_1: coordinates of one state
        :param coord_2: coordinates of another state
        :return: returns the work vector
        """
        return np.asarray(coord_1).ravel() - np.asarray(coord_2).ravel()