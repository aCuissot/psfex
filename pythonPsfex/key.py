import numpy as np
from define import * 
from enum import Enum

class cattypeenum(Enum):
    P_FLOAT = 1
    P_INT = 2
    P_STRING = 3
    P_BOOL = 4
    P_KEY = 5
    P_INTLIST = 6
    P_FLOATLIST = 7
    P_BOOLLIST = 8
    P_KEYLIST = 9
    P_STRINGLIST = 10
    
FIND_STRICT = 0
FIND_NOSTRICT = 1

class pkeystruct():
    def __init__(self, name, type, ptr, imin, imax, dmin, pc, npc, nlistmin,
                 nlistmax,nlistptrn, flag):
        self.name = name
        self.type = type
        self.ptr = ptr
        self.imin = imin
        self.imax = imax
        self.dmin = dmin
        self.dmax = dmax
        self.keylist = keylist
        self.nlistmin = nlistmin
        self.nlistmax = nlistmax
        self.nlistptr = nlistptr
        self.flag = flag

pkeydt = np.dtype([('chars', np.byte, 32*33), ('ints', np.int32, 10), ('doubles', np.float64, 2)])