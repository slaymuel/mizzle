#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
"""Wrapper for the Wetter module which performs the actual algorithm. :)

Example
-------
Usage::

    $ ./main.py TiO2_110.pdb -c config.wet

Example config file::

    atom Ti: surface
    water: 1.0
    hydroxyl: 0.0

    atom Ti: defect
        water: 0.5
        hydroxyl: 0.5
    end

Notes
-----
It is also possible to directly import the Wetter module.
 
    Example usage::

        from drizzle import Wetter

        wet = Wetter(args.verbose, newtopol)
        wet.remove_low_coordinated(Nmax, 'element')
        wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 1,\
'OH': hydroxylFrac, 'OH2': waterFrac, 'O':0.05})
        wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 2,\
'OH': fraction, 'OH2': 0, 'O':0.1})
        wet.maximize_distance()
        wet.wet()
        wet.append_atoms(fileWet, resname)
        wet.save()

"""

from radish import Topologizer
import sys
import argcomplete
import argparse
import numpy as np
import os
from IPython import embed

#"Name of program" imports
from WetParser import parse
from WetParser import get_max_coordination
from pdbExplorer import remove_low_coordinated
from pdbExplorer import append_atoms
from Wetter import Wetter
import WetParser

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def main(argv=None):

    if argv is None:
        argv = sys.argv

    argparser = argparse.ArgumentParser(description=__doc__,
                                       formatter_class=CustomFormatter)
    
    # Options
    argparser.add_argument("-s", "--silent", action="store_true",
                            help="Silence is golden")

    argparser.add_argument("coords", metavar="pdb", default="input.pdb",
                            help="Input: Structure file in pdb format")
    argparser.add_argument("-c", "--conf", default=os.path.join(os.path.dirname(WetParser.__file__), "config.wet"),
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
            print (" " + atom.get('element') + "\t    " + \
                   str(atom.get('coordination')) + "\t     " +\
                   str(hydroxylFrac) + "\t      " + str(waterFrac) +\
                   "\t         " + str(Nmax))
        elif coordination == "defect":
            dHydroxylFrac = float(atom.get('hydroxyl', None))
            dWaterFrac = float(atom.get('water', None))
            print (" " + atom.get('element') + "\t    " + \
                   str(atom.get('coordination')) + "\t     " + \
                   str(dHydroxylFrac) + "\t      " + str(dWaterFrac) + \
                   "\t         " + str(Nmax))
        else:
            raise ValueError("Unknown keyword {}".format(coordination))

    print("Solvate residues will be named: " + resname)

    ### RUN ALGORITHM ###

    # Create Wetter object
    wet = Wetter(args.coords, silent = args.silent)
    wet.remove_low_coordinated(Nmax, element)

    # Specify which types of atoms that should be hydrated
    wet.solvate({"Nmax": Nmax, "element": element, "coordination": Nmax - 1,\
                 "OH": hydroxylFrac, "OH2": waterFrac, "O": 0.05})
    wet.solvate({"Nmax": Nmax, "element": element, "coordination": Nmax - 2,\
                 "OH": dHydroxylFrac, "OH2": dWaterFrac, "O": 0.1})

    # Run minimization
    wet.optimize()

    #Create atoms and append to file
    wet.wet()
    wet.save(fileWet, resname)

if __name__ == "__main__":
    sys.exit(main())