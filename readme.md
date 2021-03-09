# Generate crystal lattices in LAMMPS or XYZ format

Writes either cubic or Al<sub>2</sub>Cu lattices to either an [xyz file](https://en.wikipedia.org/wiki/XYZ_file_format)
 or a [lammps data file](https://lammps.sandia.gov/doc/read_data.html#format-of-a-data-file).

Command line arguments:
-o: The format of the output, 'xyz' or 'lammps', default is 'xyz'
-s: The integer number of unit cells in each axis, default is 5
-d: The number density of the lattice, default = 1
-t: The lattice type, 'cubic' or 'al2cu', default is"al2cu"
-r: Ratio of A to B particles. This argument is ignored for Al2Cu lattices, default = 4

## Citation

If you use this code in work that leads to a publication please cite this code: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4592537.svg)](https://doi.org/10.5281/zenodo.4592537)

[Software is a valuable research output!](https://www.software.ac.uk/about/manifesto) 