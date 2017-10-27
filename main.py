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

from Wetter import Wetter
from radish import Topologizer
import sys
import argcomplete
import argparse
from WetParser import parse
from WetParser import get_max_coordination
from pdbExplorer import remove_lower_coordinated
from pdbExplorer import append_atoms
from shutil import copyfile
from timeit import default_timer as timer

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
    
    argparser.add_argument("Config file (default: config.wet)", nargs="+",\
                           default="config.wet", metavar="conf",\
                           help="Config file")

    argparser.add_argument("files", nargs="+", metavar="pdb",\
                           help="Input structure file in pdb format")

    # Parse arguments
    argcomplete.autocomplete(argparser)
    args = argparser.parse_args()

    # Run program
    print "Are we loud and noisy? " + str(args.verbose)
    print "The file list is:"
    print "\"" + " ".join(args.files) + "\"" 
    
    ## PREPARE INPUT ##

    #Copy pdb input file to inputfile_wet.pdb
    fileWet = args.files[0].rsplit('.', 1)[0]
    fileExt = args.files[0].rsplit('.', 1)[1]
    fileWet = fileWet + "_wet." + fileExt
    copyfile(args.files[0], fileWet)

    topol = Topologizer.from_coords(args.files[0])

    atoms = parse('config.wet')	#Call WetParser to parse config file

    element = atoms[0].get('element', None) #Get element from config.wet
    waterFrac = None
    hydroxylFrac = None
    fraction = None

    for atom in atoms:
        coordination = atom.get('coordination', None)
        if(coordination == 'high coordinated'):
            hydroxylFrac = float(atom.get('hydroxyl', None))
            waterFrac = float(atom.get('water', None))

        elif(coordination == 'low coordinated'):
            fraction = float(atom.get('fraction', None))

    Nmax = get_max_coordination(element)	#Get Nmax
    print("Nmax is: " + str(Nmax))

    #Remove reactive atoms with low coordination ( > Nmax - 2)
    newtopol = remove_lower_coordinated(topol, Nmax, element, args.verbose)

    # Save new pdb file (remove later but needed for now by pdbExplorer)
    newtopol.trj.save(fileWet, force_overwrite = True)

    ### RUN ALGORITHM ###

    wet = Wetter(args.verbose, newtopol)	#Create Wetter object
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