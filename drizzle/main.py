#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
"""Wrapper for the Wetter module which performs the actual algorithm. :)

Example
-------
Usage::

    $ ./main.py TiO2_110.pdb -c config.wet

Example config file::

    atom Ti: high coordinated
    water: 0.5
    hydroxyl: 0.5

    atom Ti: low coordinated
        fraction: 1
    end

Notes
-----
It is also possible to directly import the Wetter module.
 
    Example usage::

        from Wetter import Wetter
        from radish import Topologizer

        wet = Wetter(args.verbose, newtopol)
        wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 1,\
'OH': hydroxylFrac, 'OH2': waterFrac, 'O':0.05})
        wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 2,\
'OH': fraction, 'OH2': 0, 'O':0.1})
        wet.maximize_distance()
        wet.wet()
        we.save()

"""

from radish import Topologizer
import sys
import argcomplete
import argparse
from shutil import copyfile
from timeit import default_timer as timer
import numpy as np
import os
from IPython import embed

#"Name of program" imports
from WetParser import parse
from WetParser import get_max_coordination
from pdbExplorer import remove_low_coordinated
from pdbExplorer import append_atoms
from Wetter import Wetter

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def main(argv=None):

    if argv is None:
        argv = sys.argv

    argparser = argparse.ArgumentParser(description=__doc__,
                                       formatter_class=CustomFormatter)
    
    # Options
    argparser.add_argument("-v", "--verbose", action="store_true",
                            help="Be loud and noisy")

    argparser.add_argument("coords", metavar="pdb", default="input.pdb",
                            help="Input: Structure file in pdb format")
    argparser.add_argument("-c", "--conf", default="config.wet",
                           help="Input: Config file")
    argparser.add_argument("-o", "--out", default=None,
                           help="Output: Hydrated structure")

    # Parse arguments
    argcomplete.autocomplete(argparser)
    args = argparser.parse_args()
    
    ## PREPARE INPUT ##
    path = os.path.split(args.coords)[0]
    if(args.out):
        fileWet = os.path.join(path, args.out)
    else:
        root,ext = os.path.splitext(args.coords)
        fileWet = "{}{}{}".format(root, "_wet", ext)

    copyfile(args.coords, fileWet)

    topol = Topologizer.from_coords(args.coords)

    if topol.trj.unitcell_lengths is None:
        raise RuntimeError("No box vectors in PDB file")
    else:
        boxVectors = 10*topol.trj.unitcell_lengths.reshape(-1).astype(float)

    # Call WetParser to parse config file
    atoms, resname = parse(args.conf)

    if (not resname):
        resname = "SOL"
    # Get element from config.wet
    element = atoms[0].get("element", None)
    waterFrac = None
    hydroxylFrac = None
    fraction = None

    print("Running Wetter with Config:\n")
    print("Atom     Coordination    OH-fraction     OH2-fraction    Bulk-coordination")
    for atom in atoms:
        element = atom.get("element")
        Nmax = get_max_coordination(element)	#Get Nmax
        coordination = atom.get("coordination", None)
        if coordination == "surface":
            hydroxylFrac = float(atom.get("hydroxyl", None))
            waterFrac = float(atom.get("water", None))
        elif coordination == "defect":
            fraction = float(atom.get("fraction", None))
        else:
            raise ValueError("Unknown keyword {}".format(coordination))
        print (" " + atom.get('element') + "\t    " + \
              str(atom.get('coordination')) + "\t     " + str(hydroxylFrac) +\
               "\t      " + str(waterFrac) + "\t         " + str(Nmax))

    print("Solvate residues will be named: " + resname)
    
    # Remove reactive atoms with low coordination ( > Nmax - 2)
    newtopol = remove_low_coordinated(topol, Nmax, element, args.verbose)

    # Save new pdb file (remove later but needed for now by pdbExplorer)
    newtopol.trj.save(fileWet, force_overwrite = True)

    ### RUN ALGORITHM ###

    # Create Wetter object
    wet = Wetter(args.verbose, newtopol, boxVectors)
    # Solvate
    wet.solvate({"Nmax": Nmax, "element": element, "coordination": Nmax - 1,\
                 "OH": hydroxylFrac, "OH2": waterFrac, "O": 0.05})
    wet.solvate({"Nmax": Nmax, "element": element, "coordination": Nmax - 2,\
                 "OH": fraction, "OH2": 0, "O": 0.1})

    # Run minimization
    wet.optimize()
    
    #coords, elements = wet.wet() #Get lists of coordinates and elements
    wet.wet()
    wet.append_atoms(fileWet, resname)
    
    #Append atoms to pdbFile
    #append_atoms(file = fileWet, coords = coords, elements = elements)

if __name__ == "__main__":
    sys.exit(main())