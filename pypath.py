#!/usr/bin/env python

import argparse
from pdb.read import PDBRead
from pdb.write import PDBWrite, PDBTrajectoryWrite
from parameter.constants import Constant
from parameter.input_parameters import Parameters
from pdb.align import Align
from path.hessian import BuildHessian
from path.thermo import ThermoDynamics
from path.dynamics import Time, Transition
from report.length_check import LengthCheck
from report.fileprint import FilePrint
from report.atom_check import AtomCheck
import scipy.linalg as sp
import numpy as np
import time

start_time = time.time()

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
    file_print = FilePrint()

    ####Reading end states####

    pdb_left = PDBRead(args.start, parameters.c_alpha)
    pdb_right = PDBRead(args.end, parameters.c_alpha)

    print('Coordinates have been read\n')

    LengthCheck(pdb_left.natoms, pdb_right.natoms)

    AtomCheck(pdb_left, parameters.c_alpha)
    AtomCheck(pdb_right, parameters.c_alpha)

    ####Aligning the end states####

    align = Align(pdb_left.coord, pdb_right.coord)

    aligned_left = align.static
    aligned_right = align.moving

    flog = open('path-log', 'w')
    flog.write('Initial state = %s\n\n' % pdb_left.name)
    flog.write('Final state = %s\n\n' % pdb_right.name)
    flog.write('The RMSD between the two structures is %f\n\n' % align.rms)

    ####Building Hessian matrices####

    build_hessian = BuildHessian()

    print('@> Computing Hessian matrices.\n')

    hessian_left = build_hessian.hessian(aligned_left, pdb_left, parameters.c_alpha)
    hessian_right = build_hessian.hessian(aligned_right, pdb_right, parameters.c_alpha)

    print('\n@> The Hessian matrices have been computed.\n')

    ####Thermodynamics####

    thermo = ThermoDynamics()

    work_endpoint = thermo.work(aligned_left, aligned_right)

    energy_left = thermo.energy(hessian_left, work_endpoint)
    energy_right = thermo.energy(hessian_right, work_endpoint)

    flog.write('Energy left = %2.3f\n\n' % energy_left)
    flog.write('Energy right = %2.3f\n\n' % energy_right)

    ####Eigen decomposition####

    eval_left, evec_left = sp.eigh(hessian_left, eigvals=(0, (pdb_left.natoms * constant.dim) - 1))
    eval_right, evec_right = sp.eigh(hessian_right, eigvals=(0, (pdb_right.natoms * constant.dim) - 1))

    print('@> The Eigenvalues and the Eigenvectors have been computed\n')
    print('@> Number of modes: %d\n' % (constant.dim * pdb_left.natoms))

    if parameters.eval:
        eval_fname = args.eval

        file_print.print_array(eval_left, eval_fname + '_eval_left')
        file_print.print_array(eval_right, eval_fname + '_eval_right')

        file_print.print_multi_array(evec_left, eval_fname + '_evec_left')
        file_print.print_multi_array(evec_right, eval_fname + '_evec_right')

    path_time = Time()

    ####Transition####

    tbar_left, force_constant_left = path_time.time_to_transition_state(eval_left)
    tbar_right, force_constant_right = path_time.time_to_transition_state(eval_right)

    t_series, energy_series = path_time.time_steps(tbar_left, tbar_right, force_constant_left, force_constant_right,
                                                   energy_left, energy_right, parameters.n_conf)

    print('@> Computing transition state and trajectory.\n')

    transition = Transition(tbar_left, tbar_right, force_constant_left, force_constant_right, eval_left,
                            eval_right, evec_left, evec_right, aligned_left, aligned_right, t_series,
                            pdb_left.natoms)

    energy_left = thermo.energy(hessian_left, transition.work_left)
    energy_right = thermo.energy(hessian_right, transition.work_right)

    flog.write('The difference in energy between the two wells is %2.3f\n' % float(energy_right - energy_left))
    flog.write('\ntbar left = %.3f\n' % tbar_left)
    flog.write('\ntbar right = %.3f\n' % tbar_right)
    flog.write('\nLeft Action: %+2.3f\n' % transition.action_left)
    flog.write('\nRight Action: %+2.3f\n' % transition.action_right)
    flog.write('\nTotal Action: %+2.3f\n' % (transition.action_left + transition.action_right))

    ####Printing transition state and trajectory####

    PDBWrite(transition.xbar.reshape(pdb_left.natoms, constant.dim), pdb_left, 'trans.pdb')

    PDBTrajectoryWrite(transition.trajectory_coord.reshape(parameters.n_conf, pdb_left.natoms, constant.dim), pdb_left,
                       'trajectory.pdb')

    file_print.print_multi_array(np.column_stack((t_series, energy_series)), 'path-energy')

    print('\nTotal time taken %2.3fs\n' % (time.time() - start_time))
