import numpy as np
import Bio.PDB
import sys
from Bio.SVDSuperimposer import SVDSuperimposer
from math import *
import os

###########################
####Read a PDB file####
###########################

def pdb_read(pdb,atoms):
	if atoms == "all":
		sub_atm = []
		sub_aa = []
		sub_chain = []
		sub_aano = []
		sub_coord = []
		ca_arr = []
		j = 0
		for i in range(len(pdb)):
			if pdb[i][0:4] == 'ATOM' or pdb[i][0:6] == 'HETATM':
				sub_atm.append(pdb[i][12:16])
				sub_aa.append(pdb[i][17:20])
				sub_chain.append(pdb[i][21:22])
				sub_aano.append(pdb[i][23:26])
				sub_coord.append(float(pdb[i][30:38]))
				sub_coord.append(float(pdb[i][38:46]))
				sub_coord.append(float(pdb[i][46:54]))
				if (pdb[i][12:16].rstrip()).lstrip() == 'CA':
					ca_arr.append([])
					ca_arr[j].append(pdb[i][21:22])
					ca_arr[j].append(pdb[i][23:26])
					ca_arr[j].append(i)
					j += 1

		return sub_atm,sub_aa,sub_chain,sub_aano,sub_coord,ca_arr
	elif atoms == "CA":
		sub_atm = []
		sub_aa = []
		sub_chain = []
		sub_aano = []
		sub_coord = []
		ca_arr = []
		j = 0
		for i in range(len(pdb)):
			if pdb[i][0:4] == 'ATOM' and (pdb[i][12:16].rstrip()).lstrip() == 'CA' or pdb[i][0:6] == 'HETATM':
				sub_atm.append(pdb[i][12:16])
				sub_aa.append(pdb[i][17:20])
				sub_chain.append(pdb[i][21:22])
				sub_aano.append(pdb[i][23:26])
				sub_coord.append(float(pdb[i][30:38]))
				sub_coord.append(float(pdb[i][38:46]))
				sub_coord.append(float(pdb[i][46:54]))
				ca_arr.append([])
				ca_arr[j].append(pdb[i][21:22])
				ca_arr[j].append(pdb[i][23:26])
				ca_arr[j].append(j)
				j += 1

		return sub_atm,sub_aa,sub_chain,sub_aano,sub_coord,ca_arr

########################
####Write a PDB file####
########################

def print_pdb(sub_atm,sub_aa,sub_chain,sub_aano,sub_coord,file1):
	for i in range(len(sub_atm)):
		file1.write("ATOM   %4d %s %s %s%4d    %8.3f%8.3f%8.3f\n" % (i+1,sub_atm[i],sub_aa[i],sub_chain[i],int(sub_aano[i]),sub_coord[3*i],sub_coord[3*i+1],sub_coord[3*i+2]))
	file1.write("TER\n")
	file1.write("ENDMDL\n")
	file1.close()

########################
####Print trajectory####
########################

def print_traj_pdb(sub_atm,sub_aa,sub_chain,sub_aano,sub_coord,file1):
	for i in range(0,len(sub_coord)):
		for j in range(0,len(sub_coord[0])/3):
			file1.write("ATOM   %4d %s %s %s%4d    %8.3f%8.3f%8.3f\n" % (j+1,sub_atm[j],sub_aa[j],sub_chain[j],int(sub_aano[j]),sub_coord[i][3*j],sub_coord[i][3*j+1],sub_coord[i][3*j+2]))
		file1.write("TER\n")
		file1.write("ENDMDL\n")
	file1.close()

###############################
####Building Hessian Matrix####
###############################

def build_hessian(pdb,fconst,natm,ndim,cutoff):
	hess = [[0 for x in xrange(ndim*natm)] for x in xrange(ndim*natm)]
	for i in range(0,natm):
		for j in range(0,natm):
			if i == j:
				continue
			else:
				distsq = (pdb[3*i]-pdb[3*j])**2 + (pdb[(3*i)+1]-pdb[(3*j)+1])**2 + (pdb[(3*i)+2]-pdb[(3*j)+2])**2
				if distsq <= cutoff**2:
					for k in range(0,ndim):
						for l in range(0,ndim):
							hess[(3*i)+k][(3*j)+l] = -(fconst*(pdb[(3*i)+k]-pdb[(3*j)+k])*(pdb[(3*i)+l]-pdb[(3*j)+l]))/distsq
							hess[(3*i)+k][(3*i)+l] += -hess[(3*i)+k][(3*j)+l]
	return hess

#######################################
####Hinsen Hessian Matrix with mass####
#######################################

def build_hessian_hinsen_mass(pdb,natm,ndim,aa):
	hess = [[0 for x in xrange(ndim*natm)] for x in xrange(ndim*natm)]
	for i in range(0,natm):
		for j in range(0,natm):
			if i == j:
				continue
			else:
				distsq = (pdb[3*i]-pdb[3*j])**2 + (pdb[(3*i)+1]-pdb[(3*j)+1])**2 + (pdb[(3*i)+2]-pdb[(3*j)+2])**2
				if distsq < 16:
					fconst = 8.6e2*np.sqrt(distsq)-2.39e3
				else:
					fconst = float(128e4)/float(distsq**3)
				for k in range(0,ndim):
					for l in range(0,ndim):
						hess[(3*i)+k][(3*j)+l] = float(-(fconst*(pdb[(3*i)+k]-pdb[(3*j)+k])*(pdb[(3*i)+l]-pdb[(3*j)+l]))/distsq)/float(np.sqrt(aa_mass(aa[i]))*np.sqrt(aa_mass(aa[j])))
						hess[(3*i)+k][(3*i)+l] += -hess[(3*i)+k][(3*j)+l]
	return hess

#############################
####Constaints processing####
#############################

def constr_process(constr,ca):
	arr = []
	for i in range(0,len(constr)):
		temp = (constr[i].rstrip()).split()
		arr.append([])
		arr[i].append(float(temp[0]))
		for j in range(0,len(ca)):
			if temp[1][0:1] == ca[j][0] and temp[1][1:] == ca[j][1]:
				arr[i].append(ca[j][2])
				continue
			if temp[2][0:1] == ca[j][0] and temp[2][1:] == ca[j][1]:
				arr[i].append(ca[j][2])
				break
	return arr

##########################
####Adding Constraints####
##########################

def add_constr(hess,constr,ca,pdb,opt,cutoff,aa):
	if constr == None:
		return hess
	else:
		arr = constr_process(constr,ca)
		for idx in range(0,len(arr)):
			mat = [[ 0 for x in xrange(3)] for x in xrange(3)]
			mag = arr[idx][0]
			i = arr[idx][1]
			j = arr[idx][2]
			distsq = (pdb[3*i]-pdb[3*j])**2 + (pdb[(3*i)+1]-pdb[(3*j)+1])**2 + (pdb[(3*i)+2]-pdb[(3*j)+2])**2
			if opt  == 'CA':
				if distsq < 16:
					fconst = mag*(8.6e2*np.sqrt(distsq)-2.39e3)
				else:
					fconst = mag*(float(128e4)/float(distsq**3))
				for k in range(0,3):
					for l in range(0,3):
						hess[(3*i)+k][(3*j)+l] += float(-(fconst*(pdb[(3*i)+k]-pdb[(3*j)+k])*(pdb[(3*i)+l]-pdb[(3*j)+l]))/distsq)/float(np.sqrt(aa_mass(aa[i]))*np.sqrt(aa_mass(aa[j])))
						hess[(3*j)+k][(3*i)+l] += hess[(3*i)+k][(3*j)+l]
						hess[(3*i)+k][(3*i)+l] += -hess[(3*i)+k][(3*j)+l]
						hess[(3*j)+k][(3*j)+l] += -hess[(3*i)+k][(3*j)+l]
			elif opt == 'all':
				fconst = mag*0.01
				for k in range(0,3):
					for l in range(0,3):
						hess[(3*i)+k][(3*j)+l] += -(fconst*(pdb[(3*i)+k]-pdb[(3*j)+k])*(pdb[(3*i)+l]-pdb[(3*j)+l]))/distsq
						hess[(3*j)+k][(3*k)+l] += hess[(3*i)+k][(3*j)+l]
						hess[(3*i)+k][(3*i)+l] += -hess[(3*i)+k][(3*j)+l]
						hess[(3*j)+k][(3*j)+l] += -hess[(3*i)+k][(3*j)+l]
		return hess

##########################
####Aligning molecules####
##########################

def align(mol1,mol2,natm,ndim):
	ref_atoms = []
	alt_atoms = []
	moving2 = []

	for i in range(0,ndim*natm):
		ref_atoms.append(mol1[i])
		alt_atoms.append(mol2[i])

	fixed = np.reshape(ref_atoms,(natm,3))
	moving = np.reshape(alt_atoms,(natm,3))
	sup=SVDSuperimposer()
	sup.set(fixed,moving)
	sup.run()
	rot,tran = sup.get_rotran()
	rot=rot.astype('f')
	tran=tran.astype('f')
	for j in range(0,natm):
		moving2.append(np.dot(moving[j],rot)+tran)

	j=0

	mol2_a = [0 for x in xrange(0,natm*ndim)]

	for i in range(0,len(moving2)):
		mol2_a[j]=moving2[i][0]
		mol2_a[j+1]=moving2[i][1]
		mol2_a[j+2]=moving2[i][2]
		j = j + 3

	return ref_atoms,mol2_a,sup.get_rms()

##################
####Read files####
##################

def readfile(fname):
	fin = open(fname,"r")
	arr = fin.readlines()
	fin.close()
	return arr

###############################
####Work Matrix Calculation####
###############################

def work(array1,array2):
	sub_work = []
	for i in range(0,len(array1)):
		sub_work.append(array1[i]-array2[i])
	return sub_work

##########################
####Energy Calculation####
##########################

def ener_calc(hess,work):
	ener = 0.5*np.dot(work.T,np.dot(np.asarray(hess),work))
	return ener

####################################
####Velocity continuity equation####
####################################

def calc_xbar(eva1,eva2,evec1,evec2,mol1_d,mol2_d,tl,tr,tf,natm,ndim):
	bl = [[0 for x in xrange(natm*ndim)] for x in xrange(natm*ndim)]
	ar = [[0 for x in xrange(natm*ndim)] for x in xrange(natm*ndim)]

	for i in range(0,len(eva1)):
		if eva1[i]<0.000001:
			bl[i][i] = 1/tl
			eva2[i]=0
		else:
			if eva1[i]*tl > 707:
				bl[i][i] = eva1[i]
			else:
				bl[i][i] = (eva1[i]*cosh(eva1[i]*tl))/(sinh(eva1[i]*tl))
		if eva2[i]<0.000001:
			ar[i][i] = 1/(tl-tf)
			eva2[i]=0
		else:
			if eva2[i]*(tf-tl) > 707:
				ar[i][i] = -eva2[i]
			else:
				ar[i][i] = (eva2[i]*cosh(eva2[i]*(tl-tf)))/(sinh(eva2[i]*(tl-tf)))

	den1 = np.dot(evec1,np.dot(bl,evec1.T))
	den2 = np.dot(evec2,np.dot(ar,evec2.T))

	num1 = np.dot(den1,mol1_d)
	num2 = np.dot(den2,mol2_d)

	num = num1 - num2
	den = den1 - den2

	xbar = np.dot(num,np.linalg.inv(den))

	return xbar

########################
####Calculating tbar####
########################

def calc_tbar(eval):
	if len(eval) == 6:
		skip = 5
	else:
		skip = 6

	reci_eval = []

	for i in range(skip,len(eval)):
		reci_eval.append(float(1)/float(eval[i]))

	sum_eval = sum(reci_eval)

	incr = 0
	i = 0

	while incr <= 0.95*sum_eval:
		incr += reci_eval[i]
		i += 1

	if skip != 5:
		incr += reci_eval[i]
		i += 1

	avg_k = 1/(float(incr)/float(i))

	return (float(7)/float(avg_k)),avg_k

#################################
####Calculate path trajectory####
#################################

def path_traj(nconf,hess_mol1,hess_mol2,mol1a,mol2a,eva1,eva2,evec1,evec2,kl,kr,tl,tr,tf,workl,workr,left_e,right_e,fengy,natm,ndim):
	dt = tf/(nconf-1)
	t = 0.0

	xmat = []
	xlr = [[0 for x in xrange(natm*ndim)] for x in xrange(natm*ndim)]

	t_arr = []
	t_arr.append(0)

	h = (nconf-1)/2
	de = float(1)/float(h)

	for i in range(0,h):
		t = float(7+np.log((i+1)*de))/float(kl)
		t_arr.append(t)

	for i in range(0,h-1):
		t = float(7+np.log(1-((i+1)*de)))/float(kr)
		t_arr.append(tf-t)

	t_arr.append(tf)

	for i in range(0,nconf):
		xmat.append([])
		for j in range(0,len(eva1)):
			if t_arr[i] <= tl:
				if eva1[j]<0.000001:
					xlr[j][j] = t_arr[i]/tl
				else:
					if eva1[j]*tl > 707:
						xlr[j][j] = 1
					else:
						xlr[j][j] = (sinh(eva1[j]*t_arr[i]))/(sinh(eva1[j]*tl))
			elif t_arr[i] > tl:
				if eva2[j]<0.000001:
					xlr[j][j] = (tf-t_arr[i])/(tf-tl)
				else:
					if eva2[j]*(tf-tl) > 707:
						xlr[j][j] = -1
					else:
						xlr[j][j] = -(sinh(eva2[j]*(t_arr[i]-tf)))/(sinh(eva2[j]*(tf-tl)))
		if t_arr[i] <= tl:
			xt = np.dot(np.dot(evec1,np.dot(xlr,evec1.T)),workl) + mol1a
			fengy.write("%2.3f\t%+2.3f\n" % (t_arr[i],ener_calc(hess_mol1,np.asarray(work(xt,mol1a))) + float(right_e-left_e)))
		elif t_arr[i] > tl:
			xt = np.dot(np.dot(evec2,np.dot(xlr,evec2.T)),workr) + mol2a
			fengy.write("%2.3f\t%+2.3f\n" % (t_arr[i],ener_calc(hess_mol2,np.asarray(work(xt,mol2a)))))
		for j in range(len(xt)):
			xmat[i].append(xt[j])

	fengy.close()

	return xmat

####################################
####Calculate rocking trajectory####
####################################

def rock_traj(nconf,hess_mol,mol,eva,evec,natm,ndim,mode,exag):
	if len(eva) == 6:
		if mode == 1:
			m = 5
		else:
			print 'For a two atom system there is only 1 vibrational mode'
			sys.exit(1)
	else:
		m = 5 + mode

	kval = eva[m]

	tf = float(1)/float(kval)

	t = 0.0
	dt = float(tf)/float(nconf-1)

	xmat = []

	for i in range(0,nconf):
		xmat.append([])
		temp = [0 for x in xrange(len(mol))]
		for j in range(len(eva)):
			for k in range(m,m+1):
				temp[j] += (exag*evec[j][k]*sin(2*pi*eva[k]*t))
		for j in range(len(temp)):
			xmat[i].append(temp[j]+mol[j])
		t = t + dt

	return xmat

#######################
#####Aminoacid Mass####
#######################

def aa_mass(aa):
	mass_aa = {'ALA':71.0788,
				'CYS':103.1388,
				'ASP':115.0886,
				'GLU':129.1155,
				'PHE':147.1766,
				'GLY':57.0519,
				'HIS':137.1411,
				'ILE':113.1594,
				'LYS':128.1741,
				'LEU':113.1594,
				'MET':131.1926,
				'ASN':114.1038,
				'PRO':97.1167,
				'GLN':128.1307,
				'ARG':156.1875,
				'SER':87.0782,
				'THR':101.1051,
				'VAL':99.1326,
				'TRP':186.2132,
				'TYR':163.1760,
				'TRX':186.2132}
	for mass in mass_aa:
		if aa == mass:
			return mass_aa[mass]
