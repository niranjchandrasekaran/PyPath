# PyPath

### Program to calculate the most probable pathway connecting two equilibrium states of a protein

Given the PDB structures of two equilibrium states of a protein, the program computes the most probable path 
(minimum Onsager-Machlup action path). PyPath uses a Anisotropic Network Model based energy function to describe 
the biomolecular systems. PyPath is extremely fast and can be used to quickly determine the dynamics of protein 
and DNA molecules before venturing into more detailed molecular dynamics simulation algorithms.
 
#### Resources

The theory behind PyPath is described in detail in 
[Chandrasekaran et. al. 2017](https://aca.scitation.org/doi/10.1063/1.4976142), 
[Chandrasekaran et al. (2016)](http://scitation.aip.org/content/aca/journal/sdy/3/1/10.1063/1.4941599) 
and [Franklin et. al (2007)](http://nar.oxfordjournals.org/content/35/suppl_2/W477)

#### Dependencies

PyPath is written in python and it relies on libraries that are a part of the _NumPy_, _SciPy_ and _Biopython_ 
packages. 

The libraries can be downloaded from the following links

- [Numpy](http://www.numpy.org/)
- [Scipy](http://www.scipy.org/)
- [Biopython](http://biopython.org/wiki/Main_Page)

#### Running PyPath

PyPath can be run with the -h flag to display all the required and optional parameters

```
pypath.py -h
```
which displays the following help menu

```
optional arguments:
  -h, --help  show this help message and exit
  -start        Initial equilibrium state [PDB file] [Required]
  -end         Final equilibrium state [PDB file] [Required]
  -nconf      The number of conformations in the trajectory [default: 3]
  -calpha     If only C-alpha atoms are to be used in the simulation 
                  [default: all atom]
  -eval         print eigenvalues and eigenvectors to file

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

#### Things to remember

- The number of atoms in both the end states must be equal

- For large systems, all atom simulation is resource intensive.
