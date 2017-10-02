### Description:
Hydrates metal oxides with water, hydroxide and oxygen. 

Three main modules:
	.Wetter
	.AtomCreator
	.pdbExplorer

# Wetter
Calculates directional vectors which are the oxygen -> metal center vectors. From the sum of directional vectors of each center the coordinates where oxygen should be positioned is calculated.

# AtomCreator
The AtomCreator module holds functions for creating molecules according to specified coordinates and directional vectors. Rotational matrices are set up and aligns the molecule accordning to the direcitonal vector.


# pdbExplorer
This module handles writing and reading from the supplied .pdb file. It also removes low coordinated atoms using the 
remove_lower_coordinated function at the start of the program. 

### Installation instructions:
Install radish

### Usage:
```
python main.py config.wet pdf-file.pdb
```
optional flags:
```
-v (verbose) (do not remove lower coordinated) 
```

### Non-standard packages:
	.numpy
	.pandas
	.pyquaternion
	.radish

