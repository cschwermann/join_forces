# Join Forces
The scripts / programs in this repository are useful for the analysis of results produced by targeted molecular dynamisc simulations performed with LAMMPs (https://lammps.sandia.gov).
They were used in the analysis of the simulations in J. Keuter, C. Schwermann, A. Hepp, K. Bergander, J. Droste, M. R. Hansen, N. L. Doltsinis, C. Mück-Lichtenfeld, and F. Lips. “A highly unsaturated six-vertex amidosubstituted silicon cluster”. Chem. Sci. 11 (2020), 5895, https://dx.doi.org/10.1039/D0SC01427C.

## joinforces.f90
joinforces is an analysis tool which takes a two files containing reaction coordinates and average constraint forces as an input. The files are supposed to 
correspond to targeted molecular dynamics simulations for a transformation A-B, with the second file containing values for the backwards transformation B-A.
It interpolates both lists onto a common grid and fuses them at an "optimal" point, i.e. such that the free energy varies the least. 
For more information, see Supporting Information in the original publication.

Usage: 

    `joinforces <file1> <file2>`
   
The two arguments are the names of the files containing reaction coordinates and constraint forces.

## aver_lambda.f90
aver_lambda is a short program used to prepare the average constraint forces for joinforces. In a targeted molecular dynamics 
simulation, lammps produces an output file containing a lot of information, including reaction coordinate and constraint force
for each timestep. Here, it is assumed that the reaction coordinate is constant throughout the simulation. aver_lambda then calculates
the average constraint force for this reaction coordinate. This procedure has to be done for several simulations with fixed reaction coordinates
in order to finally calculate the free energy with joinforces.

Usage:

    `aver_lambda <file> [nstep] [end-coord]`
	
Arguments in <> are necessary, arguments in [ ] are optional.
* &lt;file> is the filename of the lammps TMD output,
* [nstep] is the optional maximum number of steps in the TMD output to be evaluated,
* [end-coord] is the optional maximum reaction coordinate. This is used to reverse the output for the back-transformation B-A, in order to have a common reaction coordinate.


Both programs can be compiled with any Fortran  compiler, e.g.

    `gfortran -o joinforces joinforces.f90`
