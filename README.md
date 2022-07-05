# HISPEC10_Geant4_Sim
Geant4 Simulations for future HISPEC10 experiments


This is Geant4 simulation exploring the feasibility of future HISPEC10 experiments at FAIR

It requires that Geant4 10.7 is installed.

The init_vis.mac macro is designed to use DAWN as the visualiser, therefore DAWN will also need to be installed.

In order to build the application create a new directory in HISPEC10_Geant4_Sim (i.e mkdir build).
Then once inside this new directory run cmake on main directory (i.e cmake PATH/HISPEC10_Geant4_Sim)

Symbolic links from the build directory named clx and dweiko will also need to be created to the executables in the CLX and DWEIKO directories.

For example:

ln -s PATH/HISPEC10_Geant4_Sim/DWEIKO PATH/HISPEC10_Geant4_Sim/build/dweiko

&

ln -s PATH/HISPEC10_Geant4_Sim/CLX PATH/HISPEC10_Geant4_Sim/build/clx

CLX and DWEIKO will need to be compiled individually with fortran 77 and fortran 90 compilers respectively. 

gfortran is reccomended as it includes them both. The makefile should detect the correct version to use in each case.

once the build directory is constructed with cmake and filled with sybolic links to the clx and dweiko executables the final installation can be made by calling make.
