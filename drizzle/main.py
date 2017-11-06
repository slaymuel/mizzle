#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
"""Wrapper for the Wetter module which performs the actual algorithm.

Example
-------
Usage::

    $ ./main.py config.wet TiO2_110.pdb

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

#'Name of program' imports
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

    argparser.add_argument("-c", "--conf",
                           help="Config file")
    argparser.add_argument("-f", "--pdb_file",
                           help="PDB file")
    # argparser.add_argument("files", nargs="+", metavar="pdb",\
    #                        help="Input structure file in pdb format")

    # Parse arguments
    argcomplete.autocomplete(argparser)
    args = argparser.parse_args()
    #verboseprint = print if verbose else lambda *a, **k: None
    # Run program
    
    ## PREPARE INPUT ##

    #Copy pdb input file to inputfile_wet.pdb
    #fileWet = args.files[0].rsplit('.', 1)[0]
    fileWet = args.pdb_file.rsplit('.', 1)[0]
    fileExt = args.pdb_file.rsplit('.', 1)[1]
    fileWet = fileWet + "_wet." + fileExt
    copyfile(args.pdb_file, fileWet)

    topol = Topologizer.from_coords(args.pdb_file)

    f = open(args.pdb_file, "r")
    content = f.readlines()
    f.close()

    for line in content:
        if('CRYST1' in line):
            tempBoxVectors = line
            boxVectors = np.array(tempBoxVectors.split()[1:4], dtype=float)
            break

    atoms = parse(args.conf)	#Call WetParser to parse config file

    element = atoms[0].get('element', None) #Get element from config.wet
    waterFrac = None
    hydroxylFrac = None
    fraction = None

    print("Running Wetter with Config:\n")
    print("Atom     Coordination    OH-fraction     OH2-fraction    Bulk-coordination")
    for atom in atoms:
        element = atom.get('element')
        Nmax = get_max_coordination(element)	#Get Nmax
        coordination = atom.get('coordination', None)
        if(coordination == 'surface'):
            hydroxylFrac = float(atom.get('hydroxyl', None))
            waterFrac = float(atom.get('water', None))

        elif(coordination == 'defect'):
            fraction = float(atom.get('fraction', None))
        print(" " + atom.get('element') + "\t    " + \
              str(atom.get('coordination')) + "\t     " + str(hydroxylFrac) +\
               "\t      " + str(waterFrac) + "\t         " + str(Nmax))


    #Remove reactive atoms with low coordination ( > Nmax - 2)
    newtopol = remove_low_coordinated(topol, Nmax, element, args.verbose)

    # Save new pdb file (remove later but needed for now by pdbExplorer)
    newtopol.trj.save(fileWet, force_overwrite = True)

    ### RUN ALGORITHM ###

    wet = Wetter(args.verbose, newtopol, boxVectors)	#Create Wetter object
    wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 1,\
                 'OH': hydroxylFrac, 'OH2': waterFrac, 'O':0.05})
    wet.solvate({'Nmax': Nmax, 'element': element, 'coordination': Nmax - 2,\
                 'OH': fraction, 'OH2': 0, 'O':0.1})
    wet.maximize_distance()	#Run minimization
    coords, elements = wet.wet() #Get lists of coordinates and elements

    #Append atoms to pdbFile
    append_atoms(file = fileWet, coords = coords, elements = elements)

if __name__ == "__main__":
    sys.exit(main())