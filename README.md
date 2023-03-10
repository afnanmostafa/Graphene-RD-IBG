# Graphene RD-IBG (Randomly Distributed Interlayer-bonded Bilayer Graphene)
### Afnan Mostafa (MS student at UMass Amherst)

repo for generating randomly-distributed interlayer-bonded bilayer graphene structures

## Files/Folders:

* example 					= 	directory for example data files generated for LAMMPS

* src 					=	source file directory that contains the following functions and scripts:

	* mainscript.m 			= 	main MATLAB script

	* generate_graphene.m 		= 	(function) generates graphene sheets (can also generate multi-layer graphene, i.e. 100 layers)

	* readtextfile.m 			= 	(function) reads a file to get spatial information of atoms

	* get_aligned_atoms.m 		= 	(function) gets aligned atoms from top and bottom layers

	* gen_rand_bot.m 			= 	(function) selects random bonds from aligned atoms

	* get_hydro.m 			= 	(function) gets coordinates of H for passivating interlayer bond-forming C atoms

	* write_equil_lmps.m 		= 	(function) writes LAMMPS input script for generating interlayer bonds

	* write_hydro_lmps.m 		= 	(function) writes the information of all passivating H

	* delete_overlap_atoms.m 	= 	(function) delets overlapping atoms

please contact/pull if you find any bugs or improvements or edge cases: afnanmostafa102@gmail.com
