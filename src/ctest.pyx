cdef extern from "test.h":

	int testOutFiles(char * progName)

def py_test() -> int:
	return testOutFiles('./psfex')
