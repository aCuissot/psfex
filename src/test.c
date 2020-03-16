#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_PLPLOT
#include PLPLOT_H
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "cplot.h"
#include "test.h"

extern const char       notokstr[];

static int setargs(char *args, char **argv)
{
	int count = 0;

	while (isspace(*args)) ++args;
	while (*args) {
		if (argv) argv[count] = args;
		while (*args && !isspace(*args)) ++args;
		if (argv && *args) *args++ = '\0';
		while (isspace(*args)) ++args;
		count++;
	}
	return count;
}

char **parsedargs(char *args, int *argc)
{
	char **argv = NULL;
	int    argn = 0;

	if (args && *args
			&& (args = strdup(args))
			&& (argn = setargs(args,NULL))
			&& (argv = malloc((argn+1) * sizeof(char *)))) {
		*argv++ = args;
		argn = setargs(args,argv);
	}

	if (args && !argv) free(args);

	*argc = argn;
	return argv;
}

void freeparsedargs(char **argv)
{
	if (argv) {
		free(argv[-1]);
		free(argv-1);
	}
}


int compareOutputFiles(FILE *fp1, FILE *fp2){
	char ch1 = getc(fp1);
	char ch2 = getc(fp2);
	int nbDiffs =0, pos = 0, line = 1;
	while (ch1 != EOF && ch2 != EOF){
		pos++;
		if (ch1 == '\n' && ch2 == '\n') {
			line++;
			pos = 0;
		}
		if (ch1 != ch2){
			nbDiff++;
			fprintf(OUTPUT, "line %d, ch1 = %c ch2 = %c\n", line, ch1, ch2)
		}
		ch1 = getc(fp1);
		ch2 = getc(fp2);
	}
	return nbDiffs;
}

int testOutFiles(char * progName){
	int argc;
	char ** argv;
	char * inParams = NULL;
	//f or the moment I put the full path of files, I've to change that
	inParams = 	" /mnt/NewHDD/PSFEx/psfex/tests/outTest.cat -c /mnt/NewHDD/PSFEx/psfex/tests/default.psfex";
	char * fullInParams = strcat(progName, inParams);
	char * path = "/mnt/NewHDD/PSFEx/psfex/tests/toCompare_psfexOut/";
	char * files_to_compare[1] = {"outTest.psf"}; // maybe test if the other output files are equals
	argv = parsedargs(inParams,&argc);

	char		**argkey, **argval,
	*str,*listbuf;
	int		a, narg, nim, ntok, opt, opt2;

	QMALLOC(argkey, char *, argc);
	QMALLOC(argval, char *, argc);

	/* Default parameters */
	prefs.command_line = argv;
	prefs.ncommand_line = argc;
	narg = nim = 0;
	listbuf = (char *)NULL;
	strcpy(prefs.prefs_name, "default.psfex");

	for (a=0; a<argc; a++) {
		if (*(argv[a]) == '-') {
			opt = (int)argv[a][1];
			if (strlen(argv[a])<4 || opt == '-') {
				opt2 = (int)tolower((int)argv[a][2]);
				if (opt == '-') {
					opt = opt2;
					opt2 = (int)tolower((int)argv[a][3]);
				}
				switch(opt) {
				case 'c':
					if (a<(argc-1)) {
						strncpy(prefs.prefs_name, argv[++a], MAXCHAR-1);
						prefs.prefs_name[MAXCHAR-1] = '\0';
					}
					break;
				case 'd':
					dumpprefs(opt2=='d' ? 1 : 0);
					exit(EXIT_SUCCESS);
					break;
				case 'v':
					printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
					exit(EXIT_SUCCESS);
					break;
				case 'h':
					exit(EXIT_SUCCESS);
					break;
				default:
					error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
				}
			}
			else {
				/*------ Config parameters */
				argkey[narg] = &argv[a][1];
				argval[narg++] = argv[++a];
			}
		}
		else {
			/*---- The input image filename(s) */
			for (; (a<argc) && (*argv[a]!='-'); a++) {
				str = (*argv[a] == '@'? listbuf=list_to_str(argv[a]+1) : argv[a]);
				for (ntok=0; (str=strtok(ntok?NULL:str, notokstr)); nim++,ntok++) {
					if (nim<MAXFILE) {
						prefs.incat_name[nim] = str;
					} else {
						error(EXIT_FAILURE, "*Error*: Too many input catalogs: ", str);
					}
				}
			}
			a--;
		}
	}

	prefs.ncat = nim;

	readprefs(prefs.prefs_name, argkey, argval, narg);
	useprefs();

	makeit();
	FILE *outputFile, *refFile;
	for (int i =0; i<1; ++i){
		outputFile = fopen(files_to_compare[i],  "r");
		char * tmpstr = strcat(path, files_to_compare[i]);
		refFile = fopen(tmpstr, "r");
		int nbDiffs = compareOutputFiles(outputFile, refFile);
		if (nbDiffs!=0){
			printf("TEST ERROR IN %s: %d bytes differs from the expected output\n", files_to_compare[i], nbDiffs);
		} else {
			printf("TEST SUCCESSFULL IN %s\n", files_to_compare[i]);

		}
	}

	exit(EXIT_SUCCESS);
}
