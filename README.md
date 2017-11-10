# Mizzle

## Motivation:
Metal oxides in the presence of water are covered in a hydration shell. To
accurately perform simulation on such systems this fact needs to be taken into
account. The motivation for this tool lies therein and thus adds hydration
shells to metal oxide systems.

## Description:
Hydrates metal oxides with OH_2, O_H and O. In `config.wet` the user specifies
which metals to hydrate and in what ratio. Hydrating a crystal containing
different metals is supported.

# Modules
The program consists of three modules:
	.Wetter
	.AtomCreator
	.pdbExplorer

## Wetter
Calculates directional vectors which are the oxygen -> metal center
vectors. From the sum of directional vectors of each center the coordinates
where oxygen should be positioned is calculated.

## AtomCreator
The AtomCreator module holds methods for creating molecules according to
specified coordinates and directional vectors from the Wetter
module. Rotational matrices are set up and aligns the molecules accordning to
the directional vectors.


## pdbExplorer
This module handles writing and reading from the supplied .pdb file. It also
handles the removal of low coordinated atoms (Nmax - 3, where Nmax is the
maximum coordination of the metal) using the remove_lower_coordinated method at
the start of the program.

# Installation

## Install using pip:

pip install .

## Development install

pip install -e .

## Rebuild Cython libraries

python setup.py build_ext --inplace


# Usage
There is one mandatory input: a file in pdb format containing the system of
interest and a config file.

```bash
mizzler input.pdb -c config.wet
```
optional flags:
```
-v (verbose) -o (ouput filename)
```
# Examples

# Other notes
## Non-standard packages:
	.numpy
	.pandas
	.pyquaternion
	.radish

