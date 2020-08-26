import numpy as np
import re
from define import *
from makeit import *
from prefs import *
from astropy.io.fits import *
from prefs import *
import sys


SYNTAX = EXECUTABLE + " catalog1 [catalog2,...][@catalog_list1 [@catalog_list2 ...]]\n" + "\t\t[-c <config_file>][-<keyword> <value>]\n" + "> to dump a default configuration file: " + EXECUTABLE + " -d \n" + "> to dump a default extended configuration file: " + EXECUTABLE + " -dd \n"

def strtok(L, R):
    return [int(x) for x in filter(None, re.split('[,' + R + ']', L))]

def main(argc, argv):
    argkey , argval = np.chararray(argc), np.chararray(argc)
    if (argc<2):
        print("\n         %s  Version %s (%s)\n" % (BANNER,MYVERSION,DATE))
        print("\nWritten by %s\n" % (AUTHORS))
        print("Copyright %s\n" % (COPYRIGHT))
        print("\nVisit %s\n" % (WEBSITE))
        print("\n%s\n" % (DISCLAIMER))
        print("Modified version 0.0.1 with test\n")
        error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX)
    
    prefs.command_line = argv
    prefs.ncommand_line = argc
    narg = nim = 0
    listbuf = None
    prefs.prefs_name = "default.psfex"

    for a in range (1 ,argc):
        if (argv[a] == '-'):
            opt = int(argv[a][1])
            if (len(argv[a])<4 or opt == '-'):
                opt2 = int(argv[a][2])
                if (opt == '-'):
                    opt = opt2
                    opt2 = int(argv[a][3])
                
                if opt == 'c':
                    if (a<(argc-1)):
                        a+=1
                        prefs.prefs_name = argv[a]
                        prefs.prefs_name[MAXCHAR-1] = '\0'
                    
                    break
                elif opt == 'd':
                    if opt2=='d':
                        dumpprefs(1)
                    else:
                        dumpprefs(0)
                    exit(EXIT_SUCCESS)
                    break
                elif opt == 'v':
                    print("%s version %s (%s)\n" % (BANNER,MYVERSION,DATE))
                    exit(EXIT_SUCCESS)
                    break
                elif opt == 'h':
                    print("\nSYNTAX: %s" % (SYNTAX))

                    exit(EXIT_SUCCESS)
                    break
                else:
                    error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX)
                
            else :
                argkey[narg] = argv[a][1]
                a=+1
                argval[narg] = argv[a]
                narg+=1

        else : 
            while((a<argc) and (argv[a]!='-')):
                a+=1
                if argv[a] == '@':
                    listbuf=argv[a]+1
                    str = listbuf
                else:
                    str = (argv[a])

                ntok=0
                while (str) :
                    nim+=1
                    ntok+=1
                    if ntok:
                        str = strtok(None, notokstr)
                    else:
                        str = strtok(str, notokstr)
                    if (nim<MAXFILE) :
                        prefs.incat_name[nim] = str
                    else:
                        error(EXIT_FAILURE, "*Error*: Too many input catalogs: ", str)
                    
            a-=1

    prefs.ncat = nim

    readprefs(prefs.prefs_name, argkey, argval, narg)
    useprefs()

    makeit()
    NFPRINTF(OUTPUT, "")
    NPRINTF(OUTPUT, "> All done (in %.1f s)\n", prefs.time_diff)

    exit(EXIT_SUCCESS)
    
if __name__ == '__main__':
    main(len(sys.argv), sys.argv)
