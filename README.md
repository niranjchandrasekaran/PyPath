# Path
---

### Program to calculate the most probable pathway connecting two equilibrium states of a protein

It also computes the trajectory along the nth normal mode of a single structure.

Download **_path.py_** and **_func.py_** to the same folder in your local computer and refer **_running-path_** for instructions on running the program

#### Dependencies

Path is written in python and it relies on a few recipes that are a part of the __NumPy__, __SciPy__ and __Biopython__ libraries. Computers that use Python 2.6 or lower might also have to install __argparse__. The version of Python on your Unix/Linux machine can be checked using the following command

```python
python -V
```
You can download the libraries, that path depends on, from the following links

- [Numpy](http://www.numpy.org/)

- [Scipy](http://www.scipy.org/)

- [Biopython](http://biopython.org/wiki/Main_Page)

- [Argparse](https://docs.python.org/3/library/argparse.html)

#### Resources

The program path is the python implementation of the PATH algorithm (previously MinActionPATH) described in [Chandrasekaran et al. (2016)](http://scitation.aip.org/content/aca/journal/sdy/3/1/10.1063/1.4941599) and [Franklin et. al (2007)](http://nar.oxfordjournals.org/content/35/suppl_2/W477)
