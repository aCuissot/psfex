import numpy as np
from define import *
from enum import Enum
from astropy.coordinates.builtin_frames.utils import norm

class vignet_operation_enum(Enum):
    VIGNET_CPY = 1
    VIGNET_ADD = 2
    VIGNET_SUB = 3
    VIGNET_MUL = 4
    VIGNET_DIV = 5

class interpolation_enum(Enum):
    INTERP_NEARESTNEIGHBOUR = 1
    INTERP_BILINEAR = 2
    INTERP_LANCZOS2 = 3
    INTERP_LANCZOS3 = 4
    INTERP_LANCZOS4 = 5

interp_kernwidth = [1,2,4,6,8]

def vignet_resample(pix1, w1, h1, pix2, w2, h2, dx, dy, step2, stepi, dgeoxpix, dgeoypix):
    pos = [0.0, 0.0]
    modnaxisn = [0, 0]
    
    if (stepi <= 0.0):
        stepi = 1.0
    dstepi = 1.0/stepi
    mx1 = w1/2.0
    mx2 = w2/2.0
    xs1 = mx1 + dx - mx2 * step2
    if (int(xs1) >= w1):
        return RETURN_ERROR
    ixs2 = 0
    if (xs1 < 0.0):
        dix2 = int(1-xs1/step2)
        if (dix2 >= w2):
            return RETURN_ERROR
        ixs2 += dix2
        xs1 += dix2 * step2;
    nx2 = int((w1 - 1 - xs1) / step2 + 1)
    ix2 = w2 - ixs2
    if (nx2>ix2):
        nx2 = ix2
    if (nx2<=0):
        return RETURN_ERROR
    my1 = h1/2.0
    my2 = h2/2.0
    ys1 = my1 + dy - my2*step2
    if (int(ys1) >= h1):
        return RETURN_ERROR
    
    iys2 = 0
    if (ys1<0.0):
        diy2 = 1-ys1/step2
        if (diy2 >= h2):
            return RETURN_ERROR
        iys2 += diy2
        ys1 += diy2*step2
    ny2 = int(((h1-1-ys1)/step2+1))
    iy2=h2-iys2
    if (ny2>iy2):
        ny2 = iy2
    if (ny2<=0):
        return RETURN_ERROR
    if (not pix2):
        pix2 = statpix2
    else:
        pix2[:w2*h2] = np.zeros(w2*h2)
        statpix2 = pix2;
    
    dgeopix_index = 0
    if (dgeoxpix and dgeoypix):
        modnaxisn[0] = w1
        modnaxisn[1] = h1
        off2 = iys2*w2 + ixs2
        pixout0 = pix2
        dgeopix0_index = off2
        dgeoxpix0 = dgeoxpix 
        dgeoypix0 = dgeoypix
        y1 = ys1 + 1.0
        for j in range(ny2):
            y1 += step2
            dgeopix_index = dgeopix0_index
            dgeoxpix0_index += w2
            x1 = xs1 + 1.0
            for i in range(nx2):
                x1 += step2
                pos[0] = x1 - dgeoxpix[dgeopix_index]*step2
                pos[1] = y1 - dgeoypix[dgeopix_index]*step2
                pixout[dgeopix_index] = vignet_interpolate_pix(pos, pix1, modnaxisn, INTERPTYPE)
                dgeopix_index += 1

        return RETURN_OK
    

    iys1a = int(ys1)
    hmh = int((INTERPW/2)/dstepi) + 2
    interph = 2*hmh
    hmw = int((INTERPW/2)/dstepi) + 2
    interpw =  2*hmw
    iys1a-= hmh
    if iys1a+hmh<0 or (iys1a< 0):
        iys1a = 0
    
    ny1 = int(ys1+ny2*step2)+interpw-hmh;   
    if (ny1>h1):
        ny1 = h1
    
    ny1 -= iys1a
    ys1 -= float(iys1a)
    mask = np.zeros(nx2*interpw, dtype=np.float64)
    nmask = np.zeros(nx2, dtype=np.int32)
    start = np.zeros(nx2, dtype=np.int32)

    x1 = xs1
    maskt = mask
    nmaskt = nmask
    startt = start
    startt_index = 0
    nmaskt_index = 0
    maskt_index = 0
    pixin_index = 0

    for j in range(nx2):
        x1+=step2
        ix1=int(x1)
        ix = ix1 - hmw
        dxm = (ix1 - x1 - hmw)*dstepi
        if (ix < 0):
            n = interpw+ix
            dxm -= ix*dstepi
            ix = 0
        else:
            n = interpw
        t=w1-ix
        if n>t:
            n=t
        startt_index = 0
        nmaskt_index = 0
        startt[startt_index] = ix
        nmaskt[nmaskt_index]= n
        startt_index+=1
        nmaskt_index += 1
        norm = 0.0
        x = dxm
        for i in range(n):
            x+=dstepi
            maskt[maskt_index]  = INTERPF(x)
            norm += maskt[maskt_index]
            maskt_index += 1
        
        if norm > 0:
            norm = 1.0/norm
        else:
            norm = dstepi
        maskt -= n
        for i in range(n):
            maskt[maskt_index]  *= norm
            maskt_index += 1
            
    pix12 = np.zeros(nx2*ny1, dtype=np.float32)
    pixin0 = pix1
    pixout0 = pix12
    pixin0_index= iys1a*w1
    pixout0_index = 0
    for k in range(ny1):
        pixin0_index+=w1
        pixout0_index+=1
        maskt = mask
        nmaskt = nmask
        nmaskt_index+=1
        startt = start
        pixout = pixout0
        pixout_index = pixout0_index
        for j in range(nx2):
            pixout_index+=ny1
            pixin = pixin0
            val = 0.0
            for i in range(nmaskt[nmaskt_index]):                
                val += maskt[maskt_index] * pixin[i+startt[j]]
                maskt_index += 1
            pixout[pixout_index] = float(val)
            nmaskt_index += 1

    y1 = ys1
    maskt = mask
    maskt_index = 0
    nmaskt = nmask
    startt = start
    for j in range(ny2):
        y1+=step2
        iy1=int(y1)
        iy = (iy1) - hmh
        dym = (iy1 - y1 - hmh)*dstepi
        if (iy < 0):
            n = interph+iy
            dym -= iy*dstepi
            iy = 0
        else:
            n = interph
        t=ny1-iy
        if (n>t):
            n=t
        startt[j] = iy
        nmaskt[j] = n
        norm = 0.0
        y=dym
        for i in range(n):
            y+=dstepi
            maskt[maskt_index] = INTERPF(y)
            norm += maskt[maskt_index]
            maskt_index += 1

        if norm > 0:
            norm = 1.0/norm
        else:
            norm = dstepi 
        maskt -= n
        for i in range (n):
            maskt[maskt_index] *= norm
            maskt_index += 1
        
    pixin0 = pix12
    pixin0_index = 0
    pixout0 = pix2
    pixout0_index = ixs2+iys2*w2
    for k in range(nx2):
        pixin0_index+=ny1
        pixout0_index+=1
        maskt = mask
        maskt_index = 0
        nmaskt = nmask
        startt = start
        pixout = pixout0
        pixout_index = 0
        for j in range(ny2):
            pixout_index+=w2
            pixin = pixin0
            val = 0.0
            for i in range(nmaskt[j]):
                val += maskt[maskt_index] * pixin[i+startt[j]]
                maskt_index += 1
            pixout[pixout_index] = float(val)
            pixout_index += 1

    return RETURN_OK


def vignet_copy(pix1, w1, h1, pix2, w2, h2, idx, idy, vigop):
    pix1_index = 0
    pix2_index = 0
    if (vigop==VIGNET_CPY):
        pix2[:w2*h2] = np.zeros(w2*h2, dtype=np.float32)
    
    ymin = h2/2+idy-h1/2
    ny=h2-ymin
    if (ny>h1):
        ny = h1
    elif ny<=0:
        return RETURN_ERROR
    
    if (ymin<0):
        pix1_index -= ymin*w1
        ny += ymin
    else :
        pix2_index += ymin*w2
    
    xmin = w2/2+idx-w1/2
    nx=w2-xmin
    if (ny>w1):
        nx = w1
    elif (nx<=0):
        return RETURN_ERROR
    
    if (xmin<0):
        pix1_index -= xmin
        nx += xmin
    else:
        pix2_index += xmin
    
    off1 = w1-nx
    off2 = w2-nx
    
    if vigop == VIGNET_CPY:
        for y in range(ny):
            pix1_index+=off1
            pix2_index+=off2
            
            for x in range(nx):
                pix2[pix2_index] = pix1[pix1_index]
                pix1_index += 1
                pix2_index += 1

    elif vigop == VIGNET_ADD:
        for y in range(ny):

            pix1_index+=off1
            pix2_index+=off2
            for x in range(nx): 
                pix2[pix2_index] += pix2[pix1_index]
                pix1_index += 1
                pix2_index += 1

    elif vigop ==  VIGNET_SUB:
        for y in range(ny):

            pix1_index+=off1
            pix2_index+=off2
            for x in range(nx):
                pix2[pix2_index] -= pix2[pix1_index]
                pix1_index += 1
                pix2_index += 1

    elif vigop ==  VIGNET_MUL:
        for y in range(ny):

            pix1_index+=off1
            pix2_index+=off2
            for x in range(nx):
                pix2[pix2_index] *= pix2[pix1_index]
                pix1_index += 1
                pix2_index += 1

    elif vigop ==  VIGNET_DIV:
        for y in range(ny):

            pix1_index+=off1
            pix2_index+=off2
            for x in range(nx, -1, -1):
                pix2_index+=1  
                if (pix1[pix1_index]):
                    pix2[pix2_index] /= pix1[pix1_index]
                    pix1_index+=1
                else:
                    if pix2[pix2_index]>0.0:
                        pix2[pix2_index] = BIG
                    else:
                        pix2[pix2_index] = -BIG

    else:
        error(EXIT_FAILURE, "*Internal Error*: unknown operation in ",
                "vignet_copy()")    

    return RETURN_OK


def vignet_aperflux(ima, var, w, h, dxc, dyc, aper, gain, backnoise, fluxvar):
    pvar = backnoise*backnoise
    if pvar>0.0:
        invbacknoise2 = 1.0/pvar
    else:
        invbacknoise2 = 0.0
        
    if gain>0.0:
        invgain = 1.0/gain
    else:
        invgain = 0.0

    raper = aper/2.0
    raper2 = raper*raper
    rintlim = raper - 0.75
    if rintlim>0.0>0.0:
        rintlim2 = rintlim**2
    else:
        rintlim2 = 0.0
    
    rextlim2 = (raper + 0.75)*(raper + 0.75)
    tv = sigtv = area = 0.0
    scaley = scalex = 1.0/APER_OVERSAMP
    scale2 = scalex*scaley
    offsetx = 0.5*(scalex-1.0)
    offsety = 0.5*(scaley-1.0)
    vthresh = BIG/2.0
    mx = dxc + float(w/2)
    my = dyc + float(h/2)

    xmin = int(mx-raper+0.499999)
    xmax = int(mx+raper+1.499999)
    ymin = int(my-raper+0.499999)
    ymax = int(my+raper+1.499999)
    if (xmin < 0 or xmax > w or ymin < 0 or ymax > h):
        return 0.0

    for y in range(ymain, ymax):
        pos = y*w + xmin
        imat_index = vart_index = pos
        imat = ima 
        vart = var
        for x in range(xmin, xmax):
            imat_index+=1
            vart_index+=1
            dx = x - mx
            dy = y - my
            r2=dx*dx+dy*dy
            if (r2 < rextlim2):
                if (r2> rintlim2):
                    dx += offsetx
                    dy += offsety
                    locarea = 0.0
                    for sy in range(APER_OVERSAMP):
                        dy+=scaley
                        dx1 = dx
                        dy2 = dy*dy
                        for sx in range(APER_OVERSAMP):
                            dx1+=scalex
                            if (dx1*dx1+dy2<raper2):
                                locarea += scale2
                            
                else:
                    locarea = 1.0
                
                area += locarea
                pvar=vart[vart_index]
                pix=imat[imat_index]
                if (pix<=-BIG or (var and pvar>=vthresh)):
                    pos = y2*w + x2
                    pix = imat[imat_index + pos]
                    x2=int(2*mx+0.49999-x)
                    y2=int(2*my+0.49999-y)
                    if (x2>=0 and x2<w and y2>=0 and y2<h and pix>-BIG):
                        if (var) :
                            pvar = var[pos]
                            if (pvar>=vthresh) :
                                pix = pvar = 0.0
                        
                    else :
                        pix = 0.0
                        if (var):
                            pvar = 0.0
                        
                tv += locarea*pix
                sigtv += locarea*pvar
                if (pix>0.0 and gain>0.0):
                    sigtv += pix*invgain*pvar*invbacknoise2
                
    if (tv>0.0):
        sigtv += tv*invgain
    

    if (fluxvar):
        fluxvar = sqrt(sigtv)
    
    return tv
    
    
def vignet_interpolate_pix(posin, pix, naxisn, interptype):
    buffer= np.zeros(INTERP_MAXKERNELWIDTH, dtype=np.float32)
    kernel= np.zeros(INTERP_MAXKERNELWIDTH, dtype=np.float32)
    dpos = [0.0, 0.0]
    kwidth = interp_kernwidth[interptype]
    start = 0
    fac = 1
    posin_ind = 0
    for n in range(2):
        val = posin
        posin_ind += 1
        width = naxisn[n]
        if interptype==INTERP_NEARESTNEIGHBOUR:
            ival = int(val-0.50001)
        else:
            ival = int(val)
        dpos[n] = val - ival;
        ival-=kwidth/2;
        if (ival<0 or ival+kwidth<=0 or ival+kwidth>width):
            return 0.0
        
        start += ival*fac
        fac *= width
    

    vignet_make_kernel(dpos[0], kernel, interptype)
    step = naxisn[0]-kwidth
    pixin = pix
    pixin_index = start
    pixout = buffer
    for j in range(kwidth):
        val = 0.0;
        kvector = kernel;
        for i in range(kwidth):
            val += kvector[i]*pixin[pixin_index]
            pixin_index+=1
        
        pixout[j] = val
        pixin_index += step;
        

    vignet_make_kernel(dpos[1], kernel, interptype)
    pixin = buffer
    val = 0.0
    kvector = kernel
    for i in range(kwidth):
        tmp_index1 += 1
        val += kvector[i]*(pixin[i])
    
    return val


def vignet_make_kernel(pos, kernel, interptype):
    kernel_index = 0
    if (interptype == INTERP_NEARESTNEIGHBOUR):
        kernel[kernel_index] = 1
    elif (interptype == INTERP_BILINEAR):
        kernel[kernel_index] = 1.0-pos
        kernel_index += 1
        kernel[kernel_index] = pos
    elif (interptype == INTERP_LANCZOS2):
        if (pos<1e-5 and pos>-1e-5):
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 1.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
        else :
            x = -PI/2.0*(pos+1.0)
            sinx1 = sinf(x)
            cosx1 = cosf(x)
            kernel[kernel_index] = sinx1/(x*x)
            val = kernel[kernel_index]
            kernel_index += 1
            x += PI/2.0
            kernel[kernel_index] = -cosx1/(x*x)
            val += kernel[kernel_index]
            x += PI/2.0
            kernel_index += 1
            kernel[kernel_index] = -sinx1/(x*x)
            val += kernel[kernel_index]
            x += PI/2.0
            kernel_index += 1
            kernel[kernel_index] = cosx1/(x*x)
            val += kernel[kernel_index]
            val = 1.0/val
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
        
    elif (interptype == INTERP_LANCZOS3):
        if (pos<1e-5 and pos>-1e-5):
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 1.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
            kernel_index += 1
            kernel[kernel_index] = 0.0
        else: 
            x = -PI/3.0*(pos+2.0)

            sinx1 = sinf(x)
            cosx1 = cosf(x)
            kernel[kernel_index] = sinx1/(x*x)
            val = kernel[kernel_index]
            kernel_index += 1

            x += PI/3.0
            sinx2=-0.5*sinx1-0.866025403785*cosx1 / (x*x)
            kernel[kernel_index] = sinx2
            val += kernel[kernel_index];
            kernel_index += 1

            x += PI/3.0
            sinx3=-0.5*sinx1+0.866025403785*cosx1/(x*x)
            kernel[kernel_index]= sinx3
            val += kernel[kernel_index] 
            kernel_index += 1

            x += PI/3.0
            kernel[kernel_index] = sinx1/(x*x)
            val += kernel[kernel_index] 
            kernel_index += 1

            x += PI/3.0
            kernel[kernel_index] = sinx2/(x*x)
            val += kernel[kernel_index]
            kernel_index += 1

            x += PI/3.0
            kernel[kernel_index]= sinx3/(x*x)
            val += kernel[kernel_index]
            val = 1.0/val
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
    elif (interptype == INTERP_LANCZOS4):
        if (pos<1e-5 and pos>-1e-5):
            kernel[kernel_index]= 0.0
            kernel_index += 1
            kernel[kernel_index]= 0.0
            kernel_index += 1
            kernel[kernel_index]= 0.0
            kernel_index += 1
            kernel[kernel_index]= 1.0
            kernel_index += 1
            kernel[kernel_index]= 0.0
            kernel_index += 1            
            kernel[kernel_index]= 0.0
            kernel_index += 1
            kernel[kernel_index]= 0.0
            kernel_index += 1
            kernel[kernel_index]= 0.0
        else:
            x = -PI/4.0*(pos+3.0)
            sinx1 = sinf(x)
            cosx1 = cosf(x)
            kernel[kernel_index] = sinx1/(x*x)
            val = kernel[kernel_index]
            kernel_index += 1
            x += PI/4.0
            sinx2=0.707106781186*(sinx1+cosx1)/(x*x)
            kernel[kernel_index] = -sinx2
            val +=kernel[kernel_index]
            kernel_index += 1
            x += PI/4.0
            kernel[kernel_index] = cosx1/(x*x)
            val += kernel[kernel_index]
            kernel_index += 1
            x += PI/4.0
            sinx3 =0.707106781186*(cosx1-sinx1)/(x*x)
            kernel[kernel_index] = -sinx3
            val += kernel[kernel_index]
            kernel_index += 1
            x += PI/4.0
            kernel[kernel_index] = -sinx1/(x*x)
            val += kernel[kernel_index];
            kernel_index += 1
            x += PI/4.0
            kernel[kernel_index] = sinx2/(x*x)
            val += kernel[kernel_index]
            kernel_index += 1
            x += PI/4.0
            kernel[kernel_index] = -cosx1/(x*x)

            val += kernel[kernel_index]
            kernel_index += 1

            x += PI/4.0
            kernel[kernel_index] = sinx3/(x*x)
            
            val += kernel[kernel_index]
            val = 1.0/val
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
            kernel_index -= 1
            kernel[kernel_index] *= val
        
    else :
        error(EXIT_FAILURE, "*Internal Error*: Unknown interpolation type in ",
                "make_kernel()")
    
    