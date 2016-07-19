# dislocMMC

## Dislocation evolution Monte Carlo simulation code for graphitic systems interfacing with the LAMMPS atomistic simulation package

dislocMMC is a code, written in C++, to simulate the creation, motion, interaction, combination and annihilation of non-basal dislocations in carbon (graphitic) systems. At the heart of the code is a unique algorithm for analysing the connectivity of an atomic configuration to identify pairs of atoms (bonds) whose rotation leads to non-basal dislocation glide and transformation. The program then performs Stone-Wales type bond rotations followed by structural minimisation (calling LAMMPS) to evolve the structure. 

<img src="http://i68.tinypic.com/30wm7nn.jpg">

The code can be run in one of two modes (specified in the setup.in file):
  1. Steepest-descent optimisation, where dislocations are moved to lower the system energy. 
  2. Monte Carlo: dislocations are evolved using the Metropolis-Hastings algorithm, with an effective temperature specified in the setup.in file. 

To run, the code requires 2 input files: setup.in and data.in. 
  setup.in: The first line states the path to run the LAMMPS optimisation. The next line then specifies: The mode (1 or 2), the maximum number of steps, the printing interval, the bond formation cutoff, the effective temperature (in K), the optimisation tolerance (in eV), and the seed for the random number generator (a positive integer). 
  data.in: The system atomic configuration in enhanced .xyz file format. The first line is the number of atoms, the second line is the orthorhombic periodic cell dimensions (xmin xmax ymin ymax zmin zmax). The remaining lines are the atomic coordinates. 

In addition, an input file is also required for the parameters of force-field and optimisation, which must be in the same directory. An example of this is in the examples directory, using the AIREBO interatomic potential. 

The code then produces three output files: 
  traj.xyz: the full atomic coordinates at each printing interval
  bonds.xyz: the coordinates of the active bond sites (dislocations) at each interval
  energy.out: the total optimised system energy at each printing interval

The initial structural analysis at the beginning of the calculation is fairly expensive, however after this the overhead of the main code is negligible compared to the LAMMPS optimisations. Therefore for large systems, the speed of the calculation can be dramatically increased using an MPI version of LAMMPS on the appropriate hardware (see the examples/ directory for an setup file.  