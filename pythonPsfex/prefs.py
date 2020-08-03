import numpy as np
from define import *

prefs = None
def dumpprefs(state):

    dp = default_prefs
    for i in range(len(dp)):
        for j in range(len(dp[0])):
        
            if (dp[i][j] != '*'):
                print("%s\n", dp[i][j])
            elif (state):
                if j != len(dp[0])-1:
                    print("%s\n",[i][j+1])
                elif i != len(dp)-1:
                    print("%s\n",[i+1][0])                    
    return




def readprefs(filename, argkey, argval, narg):
    infile = open(filename,"r")
    if (infile) == None:
        flage = 1
        warning(filename, " not found, using internal defaults")
    else:
        flage = 0

    for i in range(len(key)):
        keylist[i] = key[i].name
    
    keylist[len(key)][0] = '\0'

    argi=0
    flagc = 0
    flagd = 1
    dp = default_prefs
    warn=0
    dp_index = 0
    dpn = len(dp[0])
    line = 0
    while(True):
        if (flagd):
            if (dp[dp_index//dpn][dp_index%dpn]):
                if (dp[dp_index//dpn][dp_index%dpn]=='*'):
                    str = dp[(dp_index+1)//dpn]
                    
                else:
                    str = dp[dp_index//dpn]
                dp_index+=1
            else :
                flagd = 0
            
        
        if (not flagc and not flagd):
            if (flage or not infile.readline(line)):
                line+=1
                flagc=1

        if (flagc) :
            if (argi<narg):
                str = "%s %s" % (argkey[argi], argval[argi])
                argi+=1
            else:
                break
            
        

        keyword = strtok(str, notokstr)
        if (keyword and keyword[0]!=0 and keyword[0]!='#'):
            if (warn>=10):
                error(EXIT_FAILURE, "*Error*: No valid keyword found in ", filename)
            
            nkey = findkeys(keyword, keylist, FIND_STRICT)
            if (nkey!=RETURN_ERROR) :
                value = strtok(None, notokstr)
#ifdef    HAVE_GETENV
                dolpos=find(value, '$')
                if (value and dolpos>-1):

                    value2t = value2
                    valuet = value
                    index_valuet = 0
                    index_value2t = 0
                    
                    while (dolpos):
                        while (valuet<dolpos):
                            value2t[index_value2t] = valuet[index_valuet]
                            index_valuet += 1 
                            index_value2t +=1
                        index_valuet +=1
                        if (value[index_valuet] == '{'):
                            index_valuet+=1
                        from itertools import takewhile
                        nc = len(list(takewhile(lambda x: x not in pred, valuet[index_valuet:])))
                        envname = valuet[index_valuet:nc]
                        envname[nc] = '\0'
                        envname_index = 0
                        if (valuet[nc+index_valuet] == '}'):
                            index_valuet+=1
                        envval=getenv(envname)
                        if (not envval):
                            error(EXIT_FAILURE, "Environment variable not found: ",
                                    envname)
                        
                        while(envval[envname_index]):
                            value2t[index_value2t] = envval[envname_index]
                            envname_index += 1
                            index_value2t +=1
                        
                        while(valuet[index_valuet] and valuet[index_valuet]!='$') :
                            value2t[index_value2t] = valuet[index_valuet]
                            index_valuet+=1
                            index_value2t+=1
                        
                        if (valuet[index_valuet]):
                            dolpos = valuet
                        else:
                            dolpos = None
                            value2t[index_value2t] = '\0'

                    value = strtok(value2, notokstr)
#endif
                value_index = 0
                if (key[nkey].type == P_FLOAT):
                    if (not value or value[0]=='#'):
                        error(EXIT_FAILURE, keyword," keyword has no value!")
                    
                    if (value[value_index]=='@') :
                        value = listbuf = value[value_index+1:]
                    
                    dval = atof(value)
                    if (dval>=key[nkey].dmin and dval<=key[nkey].dmax):
                        key[nkey].ptr[0] = dval
                    else:
                        error(EXIT_FAILURE, keyword," keyword out of range")
                    
                    break

                elif key[nkey].type == P_INT:
                    if (not value or value[0]=='#'):
                        error(EXIT_FAILURE, keyword," keyword has no value!")
                    
                    if (value=='@'):
                        value[value_index] = listbuf = value[value_index+1:]
                    
                    ival = int(value)
                    if (ival>=key[nkey].imin and ival<=key[nkey].imax):
                        key[nkey].ptr[0] = ival
                    else :
                        error(EXIT_FAILURE, keyword, " keyword out of range")
                    break

                elif key[nkey].type == P_STRING:
                    if (not value or value[0]=='#'):
                        value = ""
                    
                    if (value[value_index]=='@'):
                        value = listbuf = value[value_index+1:]
                    
                    key[nkey].ptr = value
                    break
                
                elif key[nkey].type == P_BOOL:
                    if (not value or value[0]=='#'):
                        error(EXIT_FAILURE, keyword," keyword has no value!")
                    
                    if (value[value_index]=='@'):
                        value = listbuf = value[value_index+1:]
                    cp = find("yYnN", value[0])
                    if (cp>-1):
                        if cp=='y' or cp =='Y':
                            key[nkey].ptr = 1
                        else:
                            key[nkey].ptr = 0
                    else:
                        error(EXIT_FAILURE, keyword, " value must be Y or N")
                    
                    break

                elif key[nkey].type == P_KEY:
                    if (not value or value[0]=='#'):
                        error(EXIT_FAILURE, keyword," keyword has no value!")
                    
                    if (value[value_index] =='@'):
                        value = listbuf = value[value_index+1]
                    ival = findkeys(value, key[nkey].keylist,FIND_STRICT)
                    if (ival != RETURN_ERROR):
                        (key[nkey].ptr) = ival
                    else :
                        errstr = "*Error*: %s set to an unknown keyword: " % (keyword)
                        error(EXIT_FAILURE, errstr, value)
                    
                    break

                elif key[nkey].type == P_BOOLLIST:
                    if (value and value[value_index]=='@'):
                        listbuf = value[value_index+1:]
                        value = strtok(listbuf, notokstr)
                    
                    for i in range(MAXLIST):
                        if not(value and value[0]!='#'):
                            break
                        if (i>=key[nkey].nlistmax):
                            error(EXIT_FAILURE, keyword, " has too many members")
                        
                        cp = find("yYnN", value[0])
                        if (cp>-1):
                            if cp=='y' or cp =='Y':
                                key[nkey].ptr = 1
                            else:
                                key[nkey].ptr = 0
                        else:
                            error(EXIT_FAILURE, keyword, " value must be Y or N")
                    
                        value = strtok(None, notokstr)
                    
                    if (i<key[nkey].nlistmin):
                        error(EXIT_FAILURE, keyword, " list has not enough members")
                    
                    key[nkey].nlistptr= i
                    break
                elif key[nkey].type == P_INTLIST:
                    if (value and value[value_index]=='@'):
                        listbuf = value[value_index+1]
                        value = strtok(listbuf, notokstr)
                    
                    for i in range(MAXLIST):
                        if not(value and value[0]!='#'):
                            break
                        if (i>=key[nkey].nlistmax):
                            error(EXIT_FAILURE, keyword, " has too many members")
                        
                        ival = int(value)
                        if (ival>=key[nkey].imin and ival<=key[nkey].imax):
                            (key[nkey].ptr)[i] = ival
                        else :
                            error(EXIT_FAILURE, keyword, " keyword out of range")
                        
                        value = strtok(None, notokstr)
                    
                    if (i<key[nkey].nlistmin):
                        error(EXIT_FAILURE, keyword, " list has not enough members")
                    
                    key[nkey].nlistptr = i
                    break
                
                elif key[nkey].type == P_FLOATLIST:
                    if (value and value[value_index]=='@'):
                        value = strtok(listbuf = list_to_str(value+1), notokstr)
                    
                    for i in range(MAXLIST):
                        if not(value and value[0]!='#'):
                            break
                        if (i>=key[nkey].nlistmax) :
                            error(EXIT_FAILURE, keyword, " has too many members")
                        
                        dval = float(value)
                        if (dval>=key[nkey].dmin and dval<=key[nkey].dmax):
                            key[nkey].ptr[i] = dval
                        else:
                            error(EXIT_FAILURE, keyword, " keyword out of range")
                        
                        value = strtok(None, notokstr)
                    
                    if (i<key[nkey].nlistmin):
                        error(EXIT_FAILURE, keyword, " list has not enough members")
                    
                    key[nkey].nlistptr = i
                    break

                elif key[nkey].type == P_KEYLIST:
                    if (value and value[value_index]=='@') :
                        listbuf = value[value_index+1]
                        value = strtok(listbuf, notokstr)
                    
                    for i in range(MAXLIST):
                        if not(value and value[0]!='#'):
                            break
                        if (i>=key[nkey].nlistmax):
                            error(EXIT_FAILURE, keyword, " has too many members")
                        ival = findkeys(value, key[nkey].keylist, FIND_STRICT)
                        if (ival != RETURN_ERROR):
                            ((key[nkey].ptr))[i] = ival
                        else :
                            errstr = "*Error*: %s set to an unknown keyword: " % (keyword)
                            error(EXIT_FAILURE, errstr, value)
                        
                        value = strtok(None, notokstr)
                    
                    if (i<key[nkey].nlistmin):
                        error(EXIT_FAILURE, keyword, " list has not enough members")
                    
                    key[nkey].nlistptr = i
                    break

                elif key[nkey].type == P_STRINGLIST:
                    if (value and value[value_index]=='@') :
                        listbuf = value[value_index+1]
                        value = strtok(None, notokstr)
                    
                    if (not(value) or value[0]=='#'):
                        value = ""
                        flagz = 1
                    else:
                        flagz = 0
                    
                    for i in range(MAXLIST):
                        if not(value and value[0]!='#'):
                            break
                        if (i>=key[nkey].nlistmax):
                            error(EXIT_FAILURE, keyword, " has too many members")
#                        QMALLOC(((char **)key[nkey].ptr)[i], char, MAXCHAR);                        
                        key[nkey].ptr[i] = value
                        value = strtok(None, notokstr)
                    
                    if (i<key[nkey].nlistmin):
                        error(EXIT_FAILURE, keyword, " list has not enough members")
                    if flagz:
                        key[nkey].nlistptr = 0
                    else:
                        key[nkey].nlistptr = i
                        
                    break

                else:
                    error(EXIT_FAILURE, "*Internal ERROR*: Type Unknown",
                            " in readprefs()")
                    break
                
                key[nkey].flag = 1
            else :
                warning(keyword, " keyword unknown")
                warn+=1
            

    for i in range(key[i].name[0]):
        if (not key[i].flag):
            error(EXIT_FAILURE, key[i].name, " configuration keyword missing")
        
    if (not flage) :
        infile.close()
    return



def findkeys(str, keyw, mode):
    i=0
    for i in range(len(keyw)):
        if (cistrcmp(str, keyw[i], mode)):
            return i

    return RETURN_ERROR



def cistrcmp(cs, ct, mode):
    if (mode):
        for i in range(min(len(cs), len(ct))):
            diff=cs[i]-ct[i]
            if diff:
                return diff
            
    else:
        for i in range(max(len(cs), len(ct))):
            diff=cs[i] - tolowerct[i]
            if diff :
                return diff
            

    return 0

def list_to_str(listname):
    fp=open(listname,"r")
    if (not fp) :
        error(EXIT_FAILURE, "*Error*: File not found: ", listname)

    listbuf = fp.readlines()
    
    fp.close()

    return listbuf



def useprefs():
    ashort=1
    bswapflag = ashort

    if (prefs.npsf_size<2):
        prefs.psf_size[1] = prefs.psf_size[0]
    

    if (prefs.npsf_pixsize<2):
        prefs.psf_pixsize[1] = prefs.psf_pixsize[0]
    
    if (not prefs.autoselect_flag and not prefs.psf_step):
        warning("SAMPLE_AUTOSELECT set to N and PSF_SAMPLING set to 0.0:\n",
                " PSF_SAMPLING will default to 1 pixel")
    
    if ((prefs.group_deg[0]!=0) and prefs.ncontext_group != prefs.ncontext_name):
        error(EXIT_FAILURE, "*Error*: PSFVAR_GROUPS and PSFVAR_KEYS do not ",
                "match")
    
    if (not prefs.group_deg[0]):
        prefs.ncontext_group = 0
    
    for i in range(prefs.ncontext_group):
        if (prefs.context_group[i]>prefs.ngroup_deg):
            error(EXIT_FAILURE, "*Error*: PSFVAR_GROUPS out of range for ",
                    prefs.context_name[i])
    
    if (not prefs.ncontext_group):
        prefs.ngroup_deg = 0
    
    if (prefs.newbasis_type == NEWBASIS_PCAINDEPENDENT
            and prefs.psf_mef_type == PSF_MEF_COMMON):
        prefs.newbasis_type = NEWBASIS_PCACOMMON
        warning("NEWBASIS_TYPE PCA_INDEPENDENT is incompatible with "
                "MEF_TYPE COMMON:\n", " NEWBASIS_TYPE forced to PCA_COMMON")

    flag = 0
    for i in range(prefs.ncontext_name):
        if (prefs.context_name[i]!= "HIDDEN?"):
            flag = 1
            break
        
    if (flag and prefs.hidden_mef_type == HIDDEN_MEF_INDEPENDENT
            and prefs.psf_mef_type == PSF_MEF_COMMON):
        prefs.hidden_mef_type = HIDDEN_MEF_COMMON
        warning("HIDDENMEF_TYPE INDEPENDENT is incompatible with "
                "MEF_TYPE COMMON:\n", " HIDDENMEF_TYPE forced to COMMON")
        
    i=len(prefs.psf_dir)-1
    pstr=prefs.psf_dir[i]
    if (i > 0 and pstr == '/'):
        pstr = '\0'
    i=len(prefs.homokernel_dir)-1
    pstr=prefs.homokernel_dir[i]
    if (i > 0 and pstr == '/') :
        pstr = '\0'
    
    if (prefs.outcat_name == "STDOUT"):
        prefs.outcatpipe_flag = 1
    
    flag = 0
    for i in range(prefs.ncheck_type):
        if (prefs.check_type[i] != PSF_NONE):
            flag = 1
        
    if (flag and prefs.ncheck_name!=prefs.ncheck_type):
        error(EXIT_FAILURE, "*Error*: CHECKIMAGE_NAME(s) and CHECKIMAGE_TYPE(s)",
                " are not in equal number")
    return
