README

GENERAL DESCRIPTION

The kinetic energy gain on a bond, tij<ci c+j> and the bond order <ci c+j> -- that generalizes the number of electrons in a bonding state minus that number in the antibonding state for a dimer -- provides a link between the DFT and the conventional chemical notation with single, double, etc bonds.
The program calculates the decomposition of the electronic kinetic energy into single bond contributions. It works on top of the Wannier90 program (W90) and any DFT code compatible with it.

PROGRAM FUNCTIONS

Finds H_k, the Fourier transform of the tight-binding hamiltonian given by the Wannier90 code.

Diagonalises H_k at every k point in a chosen grid, getting the electronic bands of the system on this k grid.

Projects each band state onto orbital pairs or bonds, obtaining thus the contribution of each bond to each band, and vice versa.

At each k point and for each bond, sums the contributions to the energy coming from every band, weighted by the occupancy (smeared by the Fermi-Dirac distribution). This way one gets the contribution to the total kinetic energy for every bond at a given k point.

Averages the contributions per bond through the whole k point grid and stores them.

PROGRAM FLOW

Obtain the tight-binding hamiltonian from the Wannier90 code and DFT. 

Copy the files 4kinFL.make4_1.sh and 4kinFL.c into a working directory.

Read the Fermi level from the DFT code and choose a k point grid. Edit the

4kinFL.make4_1.sh

file. Select the ‘one run’ or ‘reaction path’ option and comment the other one. Write a line with the path to the hamiltonian file produced by the Wannier90 code, the Fermi level read from the DFT file and your k grid.

Execute the 4kinFL.make4_1.sh program.

