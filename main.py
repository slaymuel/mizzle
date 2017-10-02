'''
usage: python main.py TiO2_110.pdb
'''
from Wetter import Wetter
import sys
import argcomplete
import argparse
from WetParser import parse
from pdbExplorer import remove_lower_coordinated

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

	#parse('config.wet')
	#Remove reactive atoms with low coordination ( > Nmax - 2) and save in temporary fila
	file = remove_lower_coordinated(args.files[1])
	print(file)
	#Instantiate the wetter module
	wetter = Wetter(file, args.verbose)

	#Run algorithm
	wetter.wet()

if __name__ == "__main__":
    sys.exit(main())