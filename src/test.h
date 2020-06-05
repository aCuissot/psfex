#ifndef _TEST_H_
#define	_TEST_H_
#include "structs.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_PLPLOT
#include PLPLOT_H
#endif

#include "makeit.h"
#include "prefs.h"
#include "cplot.h"


char **parsedargs(char *args, int *argc);

void freeparsedargs(char **argv);

int compareOutputFiles(FILE *fp1, FILE *fp2);

int testOutFiles(char * progName);

#endif //_TEST_H_
