#ifndef _TEST_H_
#define	_TEST_H_

char **parsedargs(char *args, int *argc);

void freeparsedargs(char **argv);

int compareOutputFiles(FILE *fp1, FILE *fp2);

int testOutFiles(char * progName);

#endif //_TEST_H_
