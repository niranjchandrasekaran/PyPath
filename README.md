Given the PDB structures of two equilibrium states of a biomolecule, the program computes the most probable path 
(minimum Onsager-Machlup action path). PyPath uses a Anisotropic Network Model based energy function to describe 
the biomolecular systems. PyPath is extremely fast and can be used to quickly determine the dynamics of protein 
and DNA molecules before venturing into more detailed molecular dynamics simulation algorithms.
 
#### Resources

The theory behind PyPath is described in detail in 
[Chandrasekaran and Carter (2017)](https://aca.scitation.org/doi/10.1063/1.4976142), 
[Chandrasekaran et al. (2016)](http://aca.scitation.aip.org/content/aca/journal/sdy/3/1/10.1063/1.4941599) 
and [Franklin et. al (2007)](http://nar.oxfordjournals.org/content/35/suppl_2/W477)

#### Dependencies

PyPath is written in python and it relies on several python libraries which can be downloaded and installed using the 
following command

```
conda env create --force --file environment.yml
conda activate pypath
```

The above commands requires [Anaconda](https://www.anaconda.com/products/individual) to be installed.

#### Running PyPath

PyPath can be run with the -h flag to display all the required and optional parameters

```
pypath.py -h
```
which displays the following help menu

```
optional arguments:
  -h, --help  show this help message and exit
  -start      Initial equilibrium state [PDB file] [Required]
  -end        Final equilibrium state [PDB file] [Required]
  -nconf      The number of conformations in the trajectory [default: 3]
  -calpha     If only C-alpha atoms are to be used in the simulation 
              [default: all atom]
  -torsion    If torsion potential should be included in the all atom potential [default: all atom anm]
  -eval       print eigenvalues and eigenvectors to file

```

#### Parameters

##### Equilibrium states

The two end states are input to the program using the -start and -end parameters. The end states are PDB files which 
need to have all the fields shown in the example below

```
ATOM      1  N   GLY A 703      40.667  -1.776   7.887
```

##### Trajectory

The number of frames to be computed can be specified using the -nconf flag. It should be noted that the number 
the frames specified includes the end states.

##### Atoms to simulate

By default, all atoms in the system are included in the simulation. By using the -calpha flag, only the CA atoms 
can be simulated. This is particularly useful for large systems as results from PyPath indicate that for large 
systems, CA only simulation generates results comparable to all atom simulations.

##### Potential Energy

By default the all atom simulations use the ANM potential. When the -torsion flag is used, the torsional potential is 
also included in the all atom potential. Though including the torsional potential improves the accuracy of the all atom 
potential, for large systems, computing the torsional potential can be computationally expensive.

##### Eigenvalues and Eigenvectors

Eigenvalues and eigenvectors can be useful for variety of purposes like computing the coefficient of variation 
or calculating the motion of the molecule about particular normal modes. Using the -eval flag eigenvalues and 
eigenvectors can be printed.

#### Output

Given these input parameters, PyPath generates the four output files.

_path-log_: log file containing important output parameters

_trans.pdb_: transition state PDB file.

_trajectory.pdb_: trajectory PDB file with nconf number of frames

_path-energy_: energy of each structure along the trajectory represented as time points

_missing_sidechain_atoms.txt_: this file is generated if the atoms necessary to build the Hessian are missing

#### Things to remember

- The number of atoms in both the end states must be equal

- For large systems, all atom simulation is resource intensive.

- Since PyPath works with stable equilibrium states of biomolecules, studying the dynamics of biomolecules that are unfolded may not be possible with this program
