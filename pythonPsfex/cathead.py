from key import *
from catout import *

refoutkey = [
  pkeystruct("SOURCE_NUMBER", "Source index",
    refoutsample.detindex, H_INT, T_LONG,
    "%10d", "", "meta.number", ""),
   pkeystruct("EXTENSION", "Extension index",
    refoutsample.extindex, H_INT, T_SHORT,
    "%4d", "", "meta.number", ""),
   pkeystruct("CATALOG_NUMBER", "File index",
    refoutsample.catindex, H_INT, T_LONG,
    "%7d", "", "meta.number", ""),
   pkeystruct("FLAGS_PSF", "PSFEx rejection flags",
    refoutsample.badflag, H_INT, T_LONG,
    "%7d", "", "meta.code", ""),
   pkeystruct("VECTOR_CONTEXT", "Context vector",
    refoutsample.context, H_FLOAT, T_DOUBLE,
    "%12.6g", "", "obs.param", "",
    1, refoutsample.ncontext),
   pkeystruct("X_IMAGE", "Position along x image axis",
    refoutsample.x, H_FLOAT, T_DOUBLE,
    "%11.4f", "pixel", "pos.cartesian.x", "pix"),
   pkeystruct("Y_IMAGE", "Position along y image axis",
    refoutsample.y, H_FLOAT, T_DOUBLE,
    "%11.4f", "pixel", "pos.cartesian.y", "pix"),
   pkeystruct("DELTAX_IMAGE", "Position offset along x image axis",
    refoutsample.dx, H_FLOAT, T_FLOAT,
    "%11.4f", "pixel", "pos.cartesian.x;arith.diff", "pix"),
   pkeystruct("DELTAY_IMAGE", "Position offset along y image axis",
    refoutsample.dy, H_FLOAT, T_FLOAT,
    "%11.4f", "pixel", "pos.cartesian.y;arith.diff", "pix"),
   pkeystruct("NORM_PSF", "Source (inverse) normalization factor",
    refoutsample.norm, H_FLOAT, T_FLOAT,
    "%12.6g", "count", "phot.flux;instr.det.psf", "ct"),
   pkeystruct("FWHM_PSF", "Source (inverse) normalization factor",
    refoutsample.fwhm, H_FLOAT, T_FLOAT,
    "%8.2f", "pixel", "phys.size.diameter;instr.det.psf", "pix",),
   pkeystruct("ELLIP_PSF", "Source (inverse) normalization factor",
    refoutsample.ellip, H_FLOAT, T_FLOAT,
    "%7.5f", "", "src.ellipticity;instr.det.psf", ""),
   pkeystruct("SNR_PSF", "Signal-to-Noise ratio",
    refoutsample.snr, H_FLOAT, T_FLOAT,
    "%10.4g", "", "stat.snr;instr.det.psf", ""),
   pkeystruct("CHI2_PSF", "PSF fitting chi2/d.o.f.",
    refoutsample.chi2, H_FLOAT, T_FLOAT,
    "%12.6g", "", "stat.fit.chi2;instr.det.psf", ""),
   pkeystruct("RESI_PSF", "PSF fitting normalized residuals",
    refoutsample.modresi, H_FLOAT, T_FLOAT,
    "%12.6g", "", "stat.fit.residual;instr.det.psf", ""),
  None
  ]