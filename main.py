#!/usr/bin/env python

import argparse
from pdb.read import PDBRead
from parameter.constants import Constant
from report.length_check import LengthCheck
from parameter.input_parameters import Parameters
from pdb.align import Align
from path.hessian import BuildHessian
from path.thermo import ThermoDynamics
import scipy.linalg as sp
from path.dynamics import Time, Transition
from pdb.write import PDBWrite, PDBTrajectoryWrite

parser = argparse.ArgumentParser(
    description='PATH algorithm - Compute the most probable path connecting two equilibrium states of a biomolecule')
parser.add_argument('-start', metavar='', required=True, help='Initial equilibrium state [PDB file] [Required]')
parser.add_argument('-end', metavar='', required=True, help='Final equilibrium state [PDB file] [Required]')
parser.add_argument('-nconf', metavar='', default=3, type=int,
                    help='The number of conformations in the trajectory [default: 3]')
parser.add_argument('-calpha', action='store_true',
                    help='If only C-alpha atoms are to be used in the simulation [default: all atom]')
parser.add_argument('-eval', metavar='', help='print eigenvalues and eigenvectors to file')

args = parser.parse_args()

if __name__ == '__main__':
    parameters = Parameters(args)
    constant = Constant()

    pdb_left = PDBRead(args.start, parameters.c_alpha)
    pdb_right = PDBRead(args.end, parameters.c_alpha)

    LengthCheck(pdb_left.natoms, pdb_right.natoms)

    align = Align(pdb_left.coord, pdb_right.coord)

    aligned_left = align.static
    aligned_right = align.moving

    build_hessian = BuildHessian()

    hessian_left = build_hessian.hessian(aligned_left, pdb_left, parameters.c_alpha)
    hessian_right = build_hessian.hessian(aligned_right, pdb_right, parameters.c_alpha)

    thermo = ThermoDynamics()

    work_endpoint = thermo.work(aligned_left, aligned_right)

    energy_left = thermo.energy(hessian_left, work_endpoint)
    energy_right = thermo.energy(hessian_right, work_endpoint)

    eval_left, evec_left = sp.eigh(hessian_left, eigvals=(0, (pdb_left.natoms * constant.dim) - 1))
    eval_right, evec_right = sp.eigh(hessian_right, eigvals=(0, (pdb_right.natoms * constant.dim) - 1))

    time = Time()

    tbar_left, force_constant_left = time.time_to_transition_state(eval_left)
    tbar_right, force_constant_right = time.time_to_transition_state(eval_right)

    t_series, energy_ratio = time.time_steps(tbar_left, tbar_right, force_constant_left, force_constant_right,
                                             parameters.n_conf)

    transition = Transition(tbar_left, tbar_right, force_constant_left, force_constant_right, eval_left,
                            eval_right, evec_left, evec_right, aligned_left, aligned_right, t_series,
                            pdb_left.natoms)

    energy_left = thermo.energy(hessian_left, transition.work_left)
    energy_right = thermo.energy(hessian_right, transition.work_right)

    PDBWrite(transition.xbar.reshape(pdb_left.natoms, constant.dim), pdb_left, 'trans.pdb')

    PDBTrajectoryWrite(transition.trajectory_coord.reshape(parameters.n_conf, pdb_left.natoms, constant.dim), pdb_left,
                       'trajectory.pdb')
