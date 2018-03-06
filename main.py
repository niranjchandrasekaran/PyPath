#!/usr/bin/env python

import argparse
from pdb.read import PDBRead
from parameter.constants import Constant
from report.length_check import LengthCheck
from parameter.input_parameters import Parameters
from pdb.align import Align
from path.hessian import BuildHessian

parser = argparse.ArgumentParser(description='PATH algorithm')
subparser = parser.add_subparsers(dest='option')

pathparser = subparser.add_parser('path',
                                  help='Compute the most probable path connecting two equilibrium states \
                                  of a biomolecule')
pathparser.add_argument('-start', metavar='', required=True, help='Initial equilibrium state [PDB file] [Required]')
pathparser.add_argument('-end', metavar='', required=True, help='Final equilibrium state [PDB file] [Required]')
pathparser.add_argument('-nconf', metavar='', default=3, type=int,
                        help='The number of conformations in the trajectory [default: 3]')
pathparser.add_argument('-calpha', action='store_true',
                        help='If only C-alpha atoms are to be used in the simulation [default: all atom]')
pathparser.add_argument('-cons', metavar='',
                        help='Constraints that are to be applied between C-alpha atoms of aminoacids')

rockparser = subparser.add_parser('rock', help='Rock a macromolecule structure along one of its normal mode')
rockparser.add_argument('-structure', metavar='', required=True, help='Equilibrium structure [PDB file] [Required]')
rockparser.add_argument('-nconf', metavar='', default=3, type=int,
                        help='The number of conformations in the trajectory [default: 3]')
rockparser.add_argument('-calpha', action='store_true',
                        help='If only C-alpha atoms are to be used in the simulation [default: all atom]')
rockparser.add_argument('-mode', metavar='',
                        help='The nth mode about which the molecule is rocked')
rockparser.add_argument('-cons', metavar='',
                        help='Constraints that are to be applied between C-alpha atoms of aminoacids')
rockparser.add_argument('-exag', metavar='', default=10.0, type=float,
                        help='The motion in exaggerated while rocking the structure such that the displacement \
                        from equilibrium is perceptible [default: 10]')

args = parser.parse_args()

if __name__ == '__main__':
    if args.option == 'path':
        parameters = Parameters(args)

        pdb_left = PDBRead(args.start, parameters.c_alpha)
        pdb_right = PDBRead(args.end, parameters.c_alpha)

        LengthCheck(pdb_left.natoms, pdb_right.natoms)

        align = Align(pdb_left.coord, pdb_right.coord)

        aligned_left = align.static
        aligned_right = align.moving

        build_hessian = BuildHessian()

        hessian_left = build_hessian.hessian(aligned_left, pdb_left, parameters.c_alpha)
        hessian_right = build_hessian.hessian(aligned_right, pdb_right, parameters.c_alpha)


    elif args.option == 'smooth':
        pdb = PDBRead(args.structure)

    else:
        print(parser.format_help())