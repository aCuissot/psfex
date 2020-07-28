import numpy as np
from define import *

class contextstruct():
    def __init__(self, name, group, pcflag, ncontext, degree, ngroup, pc, npc):
        self.name = name
        self.group = group
        self.pcflag = pcflag
        self.ncontext = ncontext
        self.degree = degree
        self.ngroup = ngroup
        self.pc = pc
        self.npc = npc
        
def context_init(names, group, ndim, degree, ngroup, pcexflag):
    d2=0
    npc_index =  0
    for d in range(ndim):
        pcflag = names[d] == "HIDDEN?"
        if (pcflag):
            context.npc[npc_index]
            npc_index+=1
            if ( not pcexflag):
                context.pcflag[d] = 1
            
        if (not pcexflag or not pcflag):
            context.name[d2] = names[d].copy()
            context.group[d2] = group[d]
            groupflag[group[d]-1]+=1
            d2+=1

    context.ncontext = d2
    g2 = 0
    for g in range(ngroup):
        if (groupflag[g]):
            context.degree[g2] = degree[g]
            for d2 in range(context.ncontext):
                if (context.group[d2]==g+1):
                    context.group[d2]=g2+1
            g2+=1

    context.ngroup = g2

    return context


def context_apply(context, psf, fields, ext, catindex, ncat):
    prime=[2.0,3.0,5.0,7.0,11.0,13.0,17.0,19.0,23.0,29.0,31.0,37.0,41.0,43.0,47.0,53.0,59.0,61.0]

    ncat += catindex
    for p in range(catindex, p<ncat):
        if (ext==ALL_EXTENSIONS):
            for e in range(fields[p].next):
                fields[p].psf[e] = psf_copy(psf)
        else:
            fields[p].psf[ext] = psf_copy(psf)
    if (not context.npc):
        return

    if (len(prime) < context.ncontext):
        error(EXIT_FAILURE, "*Internal Error*: ",
                "not enough prime numbers in context_apply()")
    
    poly = psf.poly
    if (poly.ndim != context.ncontext):
        error(EXIT_FAILURE, "*Internal Error*: ",
                "poly.ndim != context.ncontext in context_apply()")
    
    context2 = context_init(context.name, context.group, context.ncontext,context.degree, context.ngroup, CONTEXT_REMOVEHIDDEN)
    ncontext2 = context2.ncontext
    poly2 = poly_init(context2.group, context2.ncontext, context2.degree,context2.ngroup)

    c2 = 0
    for c in range(context.ncontext):
        if (context.pcflag[c]):
            dpos[c] = 1.0
        else: 
            dpos[c] = prime[c2]
            c2+=1
        
    poly_func(poly, dpos)
    for c in range (context2.ncontext):
        dpos[c] = prime[c]
    
    poly_func(poly2, dpos)

    polycopyflagt = polycopyflag
    polycopyflagt_index = 0
    for n2 in range(poly2.ncoeff):
        dval = poly2.basis[n2]
        for n in range(poly.ncoeff):
            polycopyflagt[polycopyflagt_index] = (fabs(poly.basis[n] - dval) < 0.1)
            polycopyflagt_index+=1

    npix = psf.size[0]*psf.size[1]
    comp2size = npix*poly2.ncoeff
    bcoeff2 = None
    bcoeff2size = psf.nbasis*poly2.ncoeff
    npc = context.npc
    for p in range(catindex, ncat):
        c2 = 0
        for c in range(context.ncontext):
            if (context.pcflag[c]):
                dpos[c] = (context.pc[p*npc+c2]-psf.contextoffset[c])/psf.contextscale[c]
                c2+=1
            else:
                dpos[c] = 1.0
        
        poly_func(poly, dpos)

        polycopyflagt = polycopyflag
        polycopyflagt_index = 0
        for n2 in range(poly2.ncoeff):
            for n in range(poly.ncoeff):
                if (polycopyflagt[polycopyflagt_index]):
                    comp = psf.comp + n*npix
                    comp2t = comp2 + n2*npix
                    dval = poly.basis[n]
                    for i in range(npix):
                        comp2t[n2*npix +i] += dval*(comp[n*npix +i])
                    
                    if (psf.basiscoeff):
                        bcoeff2[n2] += dval*psf.basiscoeff[n]
                    
                polycopyflagt_index += 1
            
        if (ext==ALL_EXTENSIONS):
            for e in range(fields[p].next):
                fields[p].psf[e] = psf2 = psf_inherit(context2, psf)
                if (e):
                    comp2 = psf2.comp.copy()
                    if (psf.basiscoeff):
                        bcoeff2 = psf2.basiscoeff.copy()
                else:
                    psf2.comp = comp2
                    if (psf.basiscoeff):
                        psf2.basiscoeff = bcoeff2
                    
        else:
            fields[p].psf[ext] = psf2 = psf_inherit(context2, psf)
            psf2.comp = comp2
            if (psf.basiscoeff):
                psf2.basiscoeff = bcoeff2

    return



def context_end(context):
    return

