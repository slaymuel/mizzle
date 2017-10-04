'''
usage: python main.py TiO2_110.pdb
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

def get_max_coordination(element):
    foundMax = False
    f = open("MaxCoordinations.lib")
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

	argparser.add_argument("files", nargs="+", default="topol.top",
                            metavar="filename", help="List of files")

    # Parse arguments
	argcomplete.autocomplete(argparser)
	args = argparser.parse_args()

    # Run program
	print "Are we loud and noisy? " + str(args.verbose)
	print "The file list is:"
	print "\"" + " ".join(args.files) + "\"" 
    
	fileWet = args.files[1].rsplit('.', 1)[0]
	fileWet = fileWet + "_wet.pdb"
	copyfile(args.files[1], fileWet)

	topol = Topologizer.from_coords(args.files[1])

	Nmax = get_max_coordination("Ti")
	print(Nmax)
	#parse('config.wet')
	#Remove reactive atoms with low coordination ( > Nmax - 2) and save in temporary fila
	newtopol = remove_lower_coordinated(topol, Nmax)

	#Save new pdb file (remove later but needed for now by pdbExplorer)
	newtopol.trj.save(fileWet, force_overwrite = True)

	#Instantiate the wetter module
	wetter = Wetter(args.verbose, newtopol, Nmax = Nmax, centers = ['Ti'])

	#Run algorithm
	coords, elements = wetter.wet()

	#Append atoms to pdbFile
	append_atoms(file = fileWet, coords = coords, elements = elements)

if __name__ == "__main__":
    sys.exit(main())