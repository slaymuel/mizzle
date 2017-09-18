'''
usage: python main.py TiO2_110.pdb
'''
from Wetter import Wetter
import sys
import argcomplete
import argparse

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

	wetter = Wetter(args.files[0], args.verbose)

	wetter.wet()

if __name__ == "__main__":
    sys.exit(main())