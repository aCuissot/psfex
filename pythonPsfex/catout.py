import numpy as np
from define import *
from enum import Enum
import context
import key

class cattypeenum(Enum):
	CAT_NONE = 1
	CAT_ASCII_HEAD = 2
	CAT_ASCII = 3
	CAT_ASCII_VOTABLE = 4
	CAT_FITS_LDAC = 5

class outsamplestruct():
	def __init__(self, detindex, extindex, catindex, badflag, context, ncontext, norm,
				fwhm, ellip, snr, x ,y, dx, dy, chi2, modresi):
		self.detindex = detindex
		self.extindex = extindex
		self.catindex = catindex
		self.badflag = badflag
		self.context = context
		self.ncontext = ncontext
		self.norm = norm
		self.fwhm = fwhm
		self.ellip = ellip
		self.snr = snr
		self.x = x
		self.y = y
		self.dx = dx
		self.dy = dy
		self.chi2 = chi2
		self.modresi = modresi
	
outsampledt = np.dtype([('ints', np.int32, 5), ('doubles', np.float64, 2 + context.MAXCONTEXT), ('floats', np.float32, 8)])

class outcatstruct():
	def __init__(self, outsample, ascfile, objtab, objkeys, buf, ncontext):
		self.outsample = outsample
		self.ascfile = ascfile
		self.objtab = objtab
		self.objkeys = objkeys
		self.buf = buf
		self.ncontext = ncontext

#tabstruc and FILE don t have dt, yet...		
outcatdt = np.dtype([('int', np.int32), ('outsample', outsampledt), ('key', key.keydt)])

def init_outcat(filename, ncontext):
	imtabtemplate = [
			"SIMPLE  =                    T / This is a FITS file",
			"BITPIX  =                    8 / ",
			"NAXIS   =                    2 / 2D data",
			"NAXIS1  =                    1 / Number of rows",
			"NAXIS2  =                    1 / Number of columns",
			"EXTEND  =                    T / This file may contain FITS extensions",
			"END                            "]
	if (prefs.outcat_type == CAT_NONE):
		return None
	outsample = outcat.outsample
	outcat.ncontext = refoutsample.ncontext = ncontext

	outcat.objtab = objtab = new_tab("LDAC_OBJECTS")
	dptr = int(outsample - refoutsample)
	for k in range(refoutkey[k].name[0]):
		objkeys[k] = refoutkey[k]
		key = objkeys+k
		key.ptr = (key.ptr + dptr)
		if key.naxis:
			key.nbytes = t_size[key.ttype]*(key.naxisn)
		else:
			key.nbytes = t_size[key.ttype]

		add_key(key,objtab, 0)
	
	if (prefs.outcat_type == CAT_ASCII_HEAD or prefs.outcat_type == CAT_ASCII or prefs.outcat_type == CAT_ASCII_VOTABLE):
		cat = None
		if (prefs.outcatpipe_flag):
			ascfile = stdout
		ascfile = open(filename, "w+")
		if (not prefs.outcatpipe_flag and not ascfile):
			error(EXIT_FAILURE,"*Error*: cannot open ", filename)
		key = objtab.key
		if (prefs.outcat_type == CAT_ASCII_HEAD and key ):
			n=1
			for i in range(objtab.nkey):
				key=key.nextkey
				if (key.unit):
					ascfile.write( "# %3d %-22.22s %-58.58s [%s]\n" % (n, key.name,key.comment, key.unit))
				else :
					ascfile.write( "# %3d %-22.22s %-58.58s\n" % (n, key.name,key.comment))
				if key.naxis:
					n += key.naxisn
				else:
					n+=1
			
		elif (prefs.outcat_type == CAT_ASCII_VOTABLE and objtab.key):
			rfilename = find(filename, '/')
			if (rfilename==-1) :
				rfilename = filename
			else:
				rfilename+=1
			
			write_xml_header(ascfile)
			ascfile.write(" <TABLE ID=\"Output_List\" name=\"%s/out\">\n" % (rfilename))
			ascfile.write("  <DESCRIPTION>Table of detections used by %s</DESCRIPTION>\n" %	(BANNER))
			ascfile.write( "  <!-- Now comes the definition of each %s parameter -.\n" % (BANNER))
			write_vo_fields(ascfile, objtab)
			ascfile.write("   <DATA><TABLEDATA>\n")
		
	else :
		ascfile = None
		cat = new_cat(1)
		init_cat(cat)
		cat.filename[:len(filename)] = filename
		if (open_cat(cat, WRITE_ONLY) != RETURN_OK):
			error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename)
		
		save_tab(cat, cat.tab)

		asctab.headnblock = 1 + (sizeof(imtabtemplate)-1)/FBSIZE;
		asctab.headbuf = np.chararray(asctab.headnblock*FBSIZE)
		asctab.headbuf = imtabtemplate[0]
		
		buf = asctab.headbuf
		for i in range(FBSIZE*asctab.headnblock):
			if (not buf[i]):
				buf[i] = ' '
			
		

		imtab = new_tab("LDAC_IMHEAD")
		key = new_key("Field Header Card")
		key.ptr = asctab.headbuf
		asctab.headbuf = NULL
		free_tab(asctab)
		key.naxis = 2
		key.naxisn = np.zeros(key.naxis, dtype=np.int32)
		key.naxisn[0] = 80
		key.naxisn[1] = fitsfind(key.ptr, "END     ")+1
		key.htype = H_STRING
		key.ttype = T_STRING
		key.nobj = 1
		key.nbytes = key.naxisn[0]*key.naxisn[1]
		add_key(key, imtab, 0)
		save_tab(cat, imtab)
		free_tab(imtab)
		objtab.cat = cat
		init_writeobj(cat, objtab, outcat.buf)
	

	outcat.ascfile = ascfile
	outcat.objkeys = objkeys

	return outcat


def write_outcat(outcat, set):

	nc = outcat.ncontext
	outsample = outcat.outsample
	samp = set.sample
	for n in range (set.nsample):
	
		#memset(outsample, 0, sizeof(outsample)); Mettre toutes les valeurs a 0
		outsample.detindex = samp[n].detindex
		outsample.extindex = samp[n].extindex + 1
		outsample.catindex = samp[n].catindex + 1
		outsample.badflag = samp[n].badflag
		for c in range(nc):
			outsample.context[c] = samp[n].context[c]
		
		outsample.x = samp.x
		outsample.y = samp.y
		outsample.dx = samp.dx
		outsample.dy = samp.dy
		outsample.norm = samp.norm
		outsample.fwhm = samp.fwhm
		outsample.ellip = samp.ellip
		outsample.snr = samp.snr
		outsample.chi2 = samp.chi2
		outsample.modresi = samp.modresi

		if (prefs.outcat_type == CAT_ASCII_HEAD or prefs.outcat_type == CAT_ASCII):
			print_obj(outcat.ascfile, outcat.objtab)
		elif (prefs.outcat_type == CAT_ASCII_VOTABLE) :
			voprint_obj(outcat.ascfile, outcat.objtab)
		else :
			write_obj(outcat.objtab, outcat.buf)
		
	return



def end_outcat(outcat):

	objtab = outcat.objtab
	if (prefs.outcat_type == CAT_ASCII_HEAD or prefs.outcat_type == CAT_ASCII) :
		if (not prefs.outcatpipe_flag) :
			ascfile.close()
		
	elif (prefs.outcat_type == CAT_ASCII_VOTABLE):
		outcat.ascfile.write( "    </TABLEDATA></DATA>\n")
		outcat.ascfile.write( "  </TABLE>\n")
		write_xml_meta(outcat.ascfile, None);
		outcat.ascfile.write( "</RESOURCE>\n")
		outcat.ascfile.write( "</VOTABLE>\n")
	else :
		end_writeobj(objtab.cat, objtab, outcat.buf)
	
	objtab.key = None
	objtab.nkey = 0

	if (objtab.cat):
		free_cat(objtab.cat, 1)
	else :
		free_tab(objtab)
	return


def write_vo_fields(file, objtab):
	if (not objtab or not objtab.key):
		return
	key=objtab.key
	for i in range(objtab.nkey):
		key=key.nextkey
		arraysize[0] = '\0'
		if (key.naxis>1):
			for d in range(key.naxis):
				if d:
					str = "%s%d" % ("x", key.naxisn[d])
				else:
					str = "%s%d" % (" arraysize=\"", key.naxisn[d])
				arraysize += str
			
			arraysize += "\""
		
		if key.ttype == T_BYTE:
			datatype = "unsignedByte"
			break
		elif key.ttype == T_SHORT:
			datatype = "short"
			break
		elif key.ttype ==  T_LONG:	
			datatype="int"
			break
		elif key.ttype ==  T_FLOAT:	
			datatype = "float"
			break
		elif key.ttype ==  T_DOUBLE:	
			datatype = "double"
			break
		else:
			error(EXIT_FAILURE,
				"*Internal Error*: Unknown datatype in ",
				"initcat()")
		
		file.write(
				"  <FIELD name=\"%s\" ucd=\"%s\" datatype=\"%s\" unit=\"%s\"%s>\n" % 
				(key.name, key.voucd, datatype,key.vounit, arraysize))
		file.write("   <DESCRIPTION>%s</DESCRIPTION>\n" % (key.comment))
		file.write( "  </FIELD>\n")
	

	return	
