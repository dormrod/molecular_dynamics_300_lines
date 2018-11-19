# Theoretically Speaking 
## Molecular Dynamics

A simple molecular dynamics code in 300 lines!

### Overview 
Following on from the podcast on molecular simulation methods
(https://soundcloud.com/theoreticallyspeaking/molecular-simulation-and-dna-origami)
we thought it would be great to give you a practical example
of how a molecular dynamics code works.

So here is possibly the simplest implementation, a 300 line python code
to perform molecular dynamics on a single molecule.
Whilst most codes would start with a fluid, we thought it would be fun
to design your own force field for more complex molecules.

You'll find example input files in ```input_files```, but please make your own,
and tweet any good ones you make to @theorypod.

Video tutorial to come!

### Requirements

* Python 3: we recommend the anaconda distribution for first-timers!
https://www.anaconda.com/download/
* VMD: for visualising the output, 
http://www.ks.uiuc.edu/Research/vmd/

### Running

The code script is called ```md.py```. 
When you run the code the following prompt should appear:
```text
Select an input file from the list:
[0]  ./input_files/Cl2.inpt
...
```
listing all available input files. 
To select one just type the corresponding number and enter.
The simulation should then run until completion.

### Output files
There are two output files. The first is ```energies.dat```, which has the format:
```text
Timestep     Kinetic_Energy     Potential_Energy     Total_Energy
0            0.0                9.439858780668843    9.439858780668843
5            1.9091103012508244 7.526579029318008    9.435689330568835
```
in kJmol^-1. The total energy should stay near constant, 
within the limits of the integrator - if not the code is broken!
A good check is that the fluctuations should be equal about the initial value.

The second is ```coordinates.xyz```, 
which has the standard xyz file format per time step (https://en.wikipedia.org/wiki/XYZ_file_format). 
This should be readable by vmd and give a trajectory with the command ```vmd coordinates.xyz```

WARNING: if you do not rename the files between simulations they will be overwritten!
### Input Files

To start with you will need an input file. 
These can be located in the directory containing the code or any subdirectory,
the script will find them as long as they have the file extension ```.inpt```

An example file is shown below:

```text

-------------------------------------------------
Molecule Name: Cl2
-------------------------------------------------
Simulation Parameters
100        Total time (ps)
0.01        Time step (ps)
100.0        Box size (Angstrom)
0.05        Coordinate/energy write frequency (ps)
-------------------------------------------------
Atom Data
2           Number of atoms
T           Add small initial random displacement (T/F)
1.0         Maximum displacement (Angstrom)
Name    Mass    X        Y        Z
Cl1     35.0    0.0      0.0      0.0
Cl2     35.0    2.0      0.0      0.0
-------------------------------------------------
Harmonic Bond Data
1           Number of bonds
Name1    Name2    R0(Angstrom)    K(Nm^-1)
Cl1      Cl2      2.0             100.0
-------------------------------------------------
```

The important things to note are as follows:
* The initial random displacement will change between runs
* The atom names must be unique
* The code will pick up some glaring errors in the file but not all,
if it is acting weirdly it might be an incorrect name in the potential information
* The bond potential is a classical harmonic oscillator, U=0.5k(r-r0)^2,
bonds are therefore fixed and cannot break!

### Extending The Code

Please play around with the code and let us know if you do anything cool.
Some ideas to start you off are:
* Different potential types (Lennard-Jones, Morse etc.)
* Multiple molecules 
* Introducing a (Andersen) thermostat/Langevin dynamics

That's all folks.