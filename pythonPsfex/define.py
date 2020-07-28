from sys import exit, stderr

VERSION = "3.x"
DATE = "2009-03-30"

BANNER = "PSFEx"
EXECUTABLE = "psfex"
MYVERSION = VERSION
COPYRIGHT = "2011-2016 IAP/CNRS/UPMC"
DISCLAIMER = BANNER + " comes with ABSOLUTELY NO WARRANTY\n" + "You may redistribute copies of " + BANNER + "\n" + "under the terms of the GNU General Public License."
AUTHORS = "Emmanuel BERTIN <bertin@iap.fr>"
WEBSITE = "http://astromatic.net/software/psfex"
INSTITUTE = "IAP http://www.iap.fr"

OUTPUT = stderr
BIG = 1e+30
MAXCHAR = 512
MAXFILE = 32768
ALL_EXTENSIONS = (-1)

#FSEEKO = fseeko
#FTELLO = ftello

RETURN_OK = 0
RETURN_ERROR = (-1)
RETURN_FATAL_ERROR = (-2)

PI = 3.1415926535898

SEEK_SET = 0
SEEK_CUR = 1

EXIT_SUCCESS = exit(0)
EXIT_FAILURE = exit(-1)
