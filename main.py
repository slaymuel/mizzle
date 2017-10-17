#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
'''
--------------------------------------
Example:


./main.py config.wet TiO2_110.pdb
--------------------------------------

'''
from Wetter import Wetter
from radish import Topologizer
import sys
import argcomplete
import argparse
from WetParser import parse
from pdbExplorer import remove_lower_coordinated
from pdbExplorer import append_atoms
from shutil import copyfile
from timeit import default_timer as timer

def get_max_coordination(element):
    foundMax = False
    f = open("MaxCoordinations.data")
    content = f.readlines()
    f.close()
    for line in content:
        if(element in line):
            colonIndex = line.index(':')
            maxCoordination = int(line[colonIndex + 1:])
            foundMax = True
            return maxCoordination
    if(not foundMax):
        raise ValueError("Could not find maximum coordination number for \"" + element + "\" in MaxCoordinations.lib")

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
	
	argparser.add_argument("Config file (default: config.wet)", nargs="+", default="config.wet",
                            metavar="conf", help="Config file")

	argparser.add_argument("files", nargs="+", metavar="pdb", help="Input structure file in pdb format")

    # Parse arguments
	argcomplete.autocomplete(argparser)
	args = argparser.parse_args()

    # Run program
	print "Are we loud and noisy? " + str(args.verbose)
	print "The file list is:"
	print "\"" + " ".join(args.files) + "\"" 
    
	fileWet = args.files[1].rsplit('.', 1)[0]
	fileExt = args.files[1].rsplit('.', 1)[1]
	fileWet = fileWet + "_wet." + fileExt
	copyfile(args.files[1], fileWet)

	topol = Topologizer.from_coords(args.files[1])

	atoms = parse('config.wet')

	element = atoms[0].get('element', None)
	waterFrac = None
	hydroxylFrac = None
	fraction = None

	for atom in atoms:
		coordination = atom.get('coordination', None)
		if(coordination == 'high coordinated'):
			hydroxylFrac = atom.get('hydroxyl', None)
			waterFrac = atom.get('water', None)

		elif(coordination == 'low coordinated'):
			fraction = atom.get('fraction', None)

	Nmax = get_max_coordination(atoms[0]['element'])
	print("Nmax is: " + str(Nmax))

	#Remove reactive atoms with low coordination ( > Nmax - 2) and save in temporary fila
	# start = timer()
	newtopol = remove_lower_coordinated(topol, Nmax, element)

	# Save new pdb file (remove later but needed for now by pdbExplorer)
	newtopol.trj.save(fileWet, force_overwrite = True)

	# Instantiate the wetter module
	wetter = Wetter(args.verbose, newtopol, Nmax = Nmax, center = element, highWaterFrac = waterFrac, highHydroxylFrac = hydroxylFrac, lowFrac = fraction)

	# #Run algorithm
	coords, elements = wetter.wet()

	# #Append atoms to pdbFile
	append_atoms(file = fileWet, coords = coords, elements = elements)
	# end = timer()
	# print("Total runtime: " + str(end-start))
if __name__ == "__main__":
    sys.exit(main())