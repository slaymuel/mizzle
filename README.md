### Motivation:
Metal oxides in the presence of water are covered in a hydration shell. To accurately perform simulation on such systems this fact needs to be taken into account. The motivation for this tool lies therein and thus adds hydration shells to arbitrary metal-oxed systems. 

### Description:
Hydrates metal oxides with OH2 and OH. In config.wet the user specifies which metals to hydrate and in what ratio. Hydrating a crystal containing different metals is supported.

The program consists of three main modules:
	.Wetter
	.AtomCreator
	.pdbExplorer
	.potential

# Wetter
Calculates directional vectors which are the oxygen -> metal center vectors. From the sum of directional vectors of each center the coordinates where oxygen should be positioned is calculated.

# AtomCreator
The AtomCreator module holds methods for creating molecules according to specified coordinates and directional vectors from the Wetter module. Rotational matrices are set up and aligns the molecules accordning to the directional vectors.


# pdbExplorer
This module handles writing and reading from the supplied .pdb file. It also handles the removal of low coordinated atoms (Nmax - 3, where Nmax is the maximum coordination of the metal) using the remove_lower_coordinated method at the start of the program. 

# pdbExplorer
Contains the objective function to be optimized.

### Installation instructions:
Clone repository and then run 
```
pip install .
```

### Usage:
There is one mandatory input: a file in pdb format containing the system of interest and a config file.
```
./main.py pdffile.pdb -c config.wet
```
optional flags:
```
-v (verbose) -o (ouput filename)
```
### Example run

### Non-standard packages:
	.numpy
	.pandas
	.pyquaternion
	.radish

