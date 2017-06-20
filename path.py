#!/usr/bin/env python

from math import *
import sys
import os.path
# import numpy as np
import time
import argparse
# import func
# import scipy.linalg as sp

start_time = time.time()

# parser = argparse.ArgumentParser(description = 'Computes the most probable pathway connecting two equilibrium states of a protein or rock a structure along the normal modes',\
# 	epilog='For more information check the running-path file')
# parser.add_argument('-ty',metavar='',type=str,default="path",help='use "rock" for rocking the structures [default:"path"]')
# parser.add_argument('-m',metavar='',type=int,default=1,help='The nth mode along which the structure should be rocked [default:1] [optional]')
# parser.add_argument('-f1',metavar='',required=True,help='Initial pdb [required]')
# parser.add_argument('-f2',metavar='',help='Final pdb [required only for "path"]')
# parser.add_argument('-n',metavar='',type=int,default='3',help='Number of conformations [default:3] [optional]')
# parser.add_argument('-ca',metavar='',type = int,default=0,help='1 for CA only [default:0 all atom] [optional]')
# parser.add_argument('-c',metavar='',default=None,help='Constraints between CA atoms of aminoacids [optional]')
# parser.add_argument('-exag',metavar='',default=10,type=float,help='The motion is exaggerated such that the displacement from equilibrium is perceptible [default:10] [optional]')

parser = argparse.ArgumentParser(description='Computes the most probable pathway connecting two equilibrium states of a protein. It can also rock a structure along one of its normal modes',\
	epilog='For more information check the running-path file')
parser.add_argument('-mpp',action='store_true',help='Use this flag for computing the most probable path')
parser.add_argument('-rock',action='store_true',help='Use this flag for rocking the structure along nth normal mode')
parser.add_argument('-l',metavar='',required=True,help='The starting structure [PDB file]')
parser.add_argument('-r',metavar='',help='Final structure [PDB file] (Required only for -path mode)')
parser.add_argument('-nconf',metavar='',default=3,type=int,help='The number of conformations in the trajectory [default: 3]')
parser.add_argument('-calpha',action='store_true',help='If only C-alpha atoms are to be used in the simulation [default: all atom]')
parser.add_argument('-cons',metavar='',help='Constraints that are to be applied between C-alpha atoms of aminoacids')
parser.add_argument('-exag',metavar='',default=10.0,type=float,help='The motion in exaggerated while rocking the structure such that the displacement from equilibrium is perceptible [default: 10]')

args = parser.parse_args()

# TODO: Figure out why -path doesn't work

if __name__ == '__main__':
	if args.mpp:
		pass
	elif args.rock:
		pass
	else:
		sys.stdout.write('\nUse the flag -mpp (for Most Probabele Path) or -rock (for Rocking the structure)\n\n')

#############################
####Functions and Classes####
#############################

class Input(object):
	def __init__(self,args):
		self.pdb1 = func.readfile(args.f1)
		self.pdb1name = args.f1
		self.mode = args.m
		self.type = args.ty
		self.nconf = args.n
		self.exag = args.exag

		if self.type == 'path':
			if args.f2:
				self.pdb2 = func.readfile(args.f2)
				self.pdb2name = args.f2
			else:
				print '\nThe second structure id required for calculating the most probable path'
				sys.exit(1)

		if args.ca == 0:
			self.opt = 'all'
			self.space = ''
		elif args.ca == 1:
			self.opt = 'CA'
			self.space = ' CA'

		if args.c != None:
			self.c = func.readfile(args.c)
		else:
			self.c = None

def open_file(opt):
	if opt == 'log':
		fname = 'path_log'
	if opt == 'trans':
		fname = 'trans.pdb'
	if opt == 'traj':
		fname = 'traj.pdb'
	elif opt == 'engy':
		fname  = 'engy'

	if os.path.isfile(fname):
		os.remove(fname)

	out = open(fname,'w')
	return out

class Const(object):
	def __init__(self):
		self.dim = 3
		self.cutoff = 8.0
		self.k = 0.01

def len_check(arr1,arr2,space):
	if len(arr1) != len(arr2):
		print "\n@> The two structures don't have the same number of atoms. Please check!\n"
		sys.exit(1)
	else:
		natm = len(arr1)/3

	if natm > 1:
		print "\n@> There are %d%s atoms in your molecule\n" % (natm,space)
	else:
		print "\n@> There is %d atom in your molecule. PATH requires at least two atoms\n" % natm
		sys.exit(1)

	return natm

#####################
####Main Function####
#####################

# if __name__ == '__main__':
#
# 	inp = Input(args)
#
# 	ftraj = open_file('traj')
#
# 	const = Const()
#
# 	###############
# 	####Rocking####
# 	###############
#
# 	if inp.type == 'rock':
#
# 		###########################
# 		####Reading Coordinates####
# 		###########################
#
# 		atm,aa,chain,aano,mol,ca = func.pdb_read(inp.pdb1,inp.opt)
#
# 		natm = len(atm)
#
# 		#############################
# 		####Building Hessian matrices
# 		#############################
#
# 		if inp.opt == 'all':
# 			hess = func.build_hessian(mol,const.k,natm,const.dim,const.cutoff)
# 		elif inp.opt == 'CA':
# 			hess = func.build_hessian_hinsen_mass(mol,natm,const.dim,aa)
#
# 		##########################
# 		####Adding Constraints####
# 		##########################
#
# 		hess_mol = func.add_constr(hess,inp.c,ca,mol,inp.opt,const.cutoff,aa)
#
# 		print "\n@> %2.3fs - The Hessian matrix has been built.\n" %(time.time()-start_time)
#
# 		####################################
# 		####Eigenvalues and Eigenvectors####
# 		####################################
#
# 		eva, evec = sp.eigh(hess_mol,eigvals=(0,(natm*const.dim)-1))
#
# 		print "@> %2.3fs - The Eigenvalues and the Eigenvectors have been calculated.\n" % (time.time()-start_time)
#
# 		##############################
# 		####Calculating trajectory####
# 		##############################
#
# 		xmat = func.rock_traj(inp.nconf,hess_mol,mol,eva,evec,natm,const.dim,inp.mode,inp.exag)
#
# 		func.print_traj_pdb(atm,aa,chain,aano,xmat,ftraj)
#
# 		print "@> %2.3fs - The trajectory has been computed.\n" % (time.time()-start_time)
#
# 	#############################
# 	####Most probable pathway####
# 	#############################
#
# 	elif inp.type == 'path':
#
# 		flog = open_file('log')
#
# 		###########################
# 		####Reading Coordinates####
# 		###########################
#
# 		atm1,aa1,chain1,aano1,mol1,ca1 = func.pdb_read(inp.pdb1,inp.opt)
# 		atm2,aa2,chain2,aano2,mol2,ca2 = func.pdb_read(inp.pdb2,inp.opt)
#
# 		natm = len_check(mol1,mol2,inp.space)
#
# 		mol1a,mol2a,rmsd = func.align(mol1,mol2,natm,const.dim)
#
# 		flog.write('Initial State = %s\n\n' % inp.pdb1name)
# 		flog.write('Final State = %s\n\n' % inp.pdb2name)
#
# 		flog.write("The rmsd between the two structures is %f Angstroms\n\n" % rmsd)
#
# 		#############################
# 		####Building Hessian matrices
# 		#############################
#
# 		if inp.opt == 'all':
# 			hess1 = func.build_hessian(mol1a,const.k,natm,const.dim,const.cutoff)
# 			hess2 = func.build_hessian(mol2a,const.k,natm,const.dim,const.cutoff)
# 		elif inp.opt == 'CA':
# 			hess1 = func.build_hessian_hinsen_mass(mol1a,natm,const.dim,aa1)
# 			hess2 = func.build_hessian_hinsen_mass(mol2a,natm,const.dim,aa2)
#
# 		##########################
# 		####Adding Constraints####
# 		##########################
#
# 		hess_mol1 = func.add_constr(hess1,inp.c,ca1,mol1a,inp.opt,const.cutoff,aa1)
# 		hess_mol2 = func.add_constr(hess2,inp.c,ca2,mol2a,inp.opt,const.cutoff,aa2)
#
# 		print "@> %2.3fs - The Hessian matrices have been built.\n" %(time.time()-start_time)
#
# 		####################################
# 		####Eigenvalues and Eigenvectors####
# 		####################################
#
# 		eva1, evec1 = sp.eigh(hess_mol1,eigvals=(0,(natm*const.dim)-1))
# 		eva2, evec2 = sp.eigh(hess_mol2,eigvals=(0,(natm*const.dim)-1))
#
# 		print "@> %2.3fs - The Eigenvalues and the Eigenvectors have been calculated.\n" % (time.time()-start_time)
#
# 		########################
# 		####Calculating tbar####
# 		########################
#
# 		tl,kl = func.calc_tbar(eva1)
# 		tr,kr = func.calc_tbar(eva2)
#
# 		tf = tl + tr
#
# 		########################
# 		####Identifying xbar####
# 		########################
#
# 		xbar = func.calc_xbar(eva1,eva2,evec1,evec2,mol1a,mol2a,tl,tr,tf,natm,const.dim)
#
# 		ftrans = open_file('trans')
# 		func.print_pdb(atm1,aa1,chain1,aano1,xbar,ftrans)
#
# 		workl = np.asarray(func.work(xbar,mol1a))
# 		workr = np.asarray(func.work(xbar,mol2a))
#
# 		###########################################
# 		####Calculating transition state energy####
# 		###########################################
#
# 		left_e = func.ener_calc(hess_mol1,workl)
# 		right_e = func.ener_calc(hess_mol2,workr)
#
# 		flog.write("The difference in energy between the two wells is : %2.3f\n" % (float(right_e-left_e)))
# 		flog.write("\ntbar left  = %.3e\n" % tl)
# 		flog.write("\ntbar right = %.3e\n" % tr)
# 		flog.write("\nUtrans  = %2.3f\n" % right_e)
#
# 		flog.close()
#
# 		##############################
# 		####Calculating trajectory####
# 		##############################
#
# 		fengy = open_file('engy')
#
# 		xmat = func.path_traj(inp.nconf,hess_mol1,hess_mol2,mol1a,mol2a,eva1,eva2,evec1,evec2,kl,kr,tl,tr,tf,workl,workr,left_e,right_e,fengy,natm,const.dim)
#
# 		func.print_traj_pdb(atm1,aa1,chain1,aano1,xmat,ftraj)
#
# 		print "@> %2.3fs - The trajectory has been computed.\n" % (time.time()-start_time)
#
# 	#################################
# 	####Time taken by the program####
# 	#################################
#
# 	print "Total time taken %2.3fs\n" % (time.time()-start_time)
