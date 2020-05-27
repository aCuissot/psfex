#ifndef _STRUCT__H_
#define _STRUCT__H_
#include "define.h"
#include	"fits/fitscat.h"
#include "wcs/poly.h"
#include <stdio.h>


typedef enum {PSF_NONE, PSF_BASIS, PSF_CHI, PSF_PROTO, PSF_RESIDUALS,
	PSF_RESIDUALS_GRID, PSF_SAMPLES, PSF_SAMPLES_GRID,
	PSF_SNAPSHOTS, PSF_SNAPSHOTS_IMRES,
	PSF_WEIGHTS, PSF_MOFFAT,PSF_SUBMOFFAT,PSF_SUBSYM, PSF_GREAT}
checkenum;

typedef enum {CAT_NONE, CAT_ASCII_HEAD, CAT_ASCII,
	CAT_ASCII_VOTABLE, CAT_FITS_LDAC}
cattypenum;

typedef enum {CPLOT_NONE, CPLOT_FWHM, CPLOT_ELLIPTICITY, CPLOT_MOFFATRESI,
	CPLOT_ASYMRESI, CPLOT_COUNTS, CPLOT_COUNTFRAC, CPLOT_CHI2,
	CPLOT_MODRESI, CPLOT_SNRVSFWHM}
cplotenum;

typedef enum {CPLOT_NULL, CPLOT_XWIN, CPLOT_TK, CPLOT_XTERM, CPLOT_PLMETA,
	CPLOT_PS, CPLOT_PSC, CPLOT_XFIG, CPLOT_LJIIP, CPLOT_LJHPGL, CPLOT_IMP,
	CPLOT_PBM, CPLOT_PNG, CPLOT_JPEG, CPLOT_PSTEX, CPLOT_AQT, CPLOT_PDF,
	CPLOT_SVG}
cplotdevenum;

typedef  enum {CELSYS_NATIVE, CELSYS_PIXEL, CELSYS_EQUATORIAL, CELSYS_GALACTIC,
	CELSYS_ECLIPTIC, CELSYS_SUPERGALACTIC}
celsysenum;

typedef enum {BASIS_NONE, BASIS_PIXEL, BASIS_GAUSS_LAGUERRE, BASIS_FILE,
		BASIS_PIXEL_AUTO}
basistypenum;

typedef  enum {VIGNET_CPY, VIGNET_ADD, VIGNET_SUB, VIGNET_MUL, VIGNET_DIV}
vigopenum;

typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}
interpenum;

typedef struct {
	cplotdevenum device;
	char *devname;
	char *extension;}
devicestruct;

typedef struct
{
	char		name[32];
	enum  {P_FLOAT, P_INT, P_STRING, P_BOOL, P_KEY, P_INTLIST, P_FLOATLIST,
		P_BOOLLIST, P_KEYLIST, P_STRINGLIST} type;
		void		*ptr;			/* Pointer to the keyword value */
		int		imin, imax;		/* Range for int's */
		double	dmin, dmax;		/* Range for floats */
		char		keylist[32][32];	/* List of keywords */
		int           nlistmin;		/* Minimum number of list members */
		int           nlistmax; 		/* Maximum number of list members */
		int		*nlistptr;		/* Ptr to store the nb of read params*/
		int		flag;
}	pkeystruct;

typedef struct sample
{
	int		detindex;		/* Detection index */
	int		extindex;		/* Extension index */
	int		catindex;		/* Catalogue index */
	float		*vig;			/* Vignette array */
	float		*vigresi;		/* Residual-map of the PSF-residuals */
	float		*vigchi;		/* Chi-map of the PSF-residuals */
	float		*vigweight;		/* Vignette-weight array */
	float		*vigdgeox;		/* Vignette diff geo x-map array */
	float		*vigdgeoy;		/* Vignette diff geo y-map array */
	float		norm;			/* Normalisation */
	float		fwhm;			/* Source FWHM (derived from ) */
	float		ellip;			/* Source ellipticity */
	float		snr;			/* Source Signal-to-Noise ratio */
	double	x,y;			/* x,y position estimate in frame */
	float		dx,dy;			/* x,y shift / vignet center */
	float		backnoise2;		/* Variance of the background noise */
	float		gain;			/* conversion factor (e-/ADU) */
	float		chi2;			/* Chi2 of the fit */
	float		modresi;		/* Residual index */
	double	*context;		/* Context vector */
	int		badflag;		/* True if sample is discarded */
}	samplestruct;

typedef struct set
{
	char		*head;			/* Table structure */
	struct sample	*sample;		/* Array of samples */
	int		nsample;		/* Number of samples in stack */
	int		ngood;			/* Number of "good" samples */
	int		nsamplemax;		/* Max number of samples in stack */
	int		*vigsize;		/* Dimensions of vignette frames */
	int		vigdim;			/* Dimensionality of the vignette */
	int		nvig;			/* Number of pixels of the vignette */
	int		ncontext;		/* Number of contexts */
	char		**contextname;		/* List of context keywords used */
	double	*contextoffset;		/* Offset to apply to context data */
	double	*contextscale;		/* Scaling to apply to context data */
	float		fwhm;			/* FWHM of the PSF core */
	int		badflags;		/* # discarded with bad SEx flags */
	int		badwflags;		/* # discarded with bad SEx flags */
	int		badimaflags;		/* # discarded with bad SEx flags */
	int		badsn;			/* # discarded with S/N too low*/
	int		badfwhmmin;		/* # discarded with size too small */
	int		badfwhmmax;		/* # discarded with size too large */
	int		badellip;		/* # discarded with too much ellipt. */
	int		badpix;			/* # discarded with too many bad pix. */
}	setstruct;

typedef struct moffat
{
	double	context[POLY_MAXDIM];	/* Context coordinates */
	float		amplitude;	/* Central amplitude */
	float		xc[2];		/* Center coordinates */
	float		fwhm_min;	/* FWHM along the minor axis */
	float		fwhm_max;	/* FWHM along the major axis */
	float		theta;		/* Position angle of the major axis / NAXIS1 */
	float		beta;		/* Moffat beta parameter */
	float		residuals;	/* Normalized residuals */
	float		symresiduals;	/* Normalized symmetry residuals */
	float		noiseqarea;	/* Noise equivalent area (pixels^2) */
	int		nsubpix;	/* Number of supersampled pixels */
}	moffatstruct;

typedef struct psf
{
	int		dim;		/* Dimensionality of the tabulated data */
	int		*size;		/* PSF dimensions */
	int		npix;		/* Total number of involved PSF pixels */
	float		*comp; 		/* Complete pix. data (PSF components) */
	float		*loc;		/* Local PSF */
	float		*resi;		/* Map of residuals */
	char		**contextname;	/* Array of context key-names */
	double	*contextoffset;	/* Offset to apply to context data */
	double	*contextscale;	/* Scaling to apply to context data */
	int		cx,cy;		/* Indices of X and Y mapping contexts */
	struct poly	*poly;		/* Polynom describing the PSF variations */
	float		pixstep;	/* Mask oversampling (pixel). */
	float		pixsize[2];	/* Effective pixel size on each axis (pixel) */
	int		samples_total;	/* Total number of detections */
	int		samples_loaded;	/* Number of detections loaded */
	int		samples_accepted;/* Number of detections accepted */
	double	chi2;		/* chi2/d.o.f. */
	float		fwhm;		/* Initial guess of the FWHM */
	int		*pixmask;	/* Pixel mask for local bases */
	float		*basis;		/* Basis vectors */
	float		*basiscoeff;	/* Basis vector coefficients */
	int		nbasis;		/* Number of basis vectors */
	int		ndata;		/* Size of the design matrix along data axis */
	int		nsnap;		/* Total number of snapshots */
	int		nmed;		/* Median position amongst snapshots */
	int		nsubpix;	/* Number of intrapixel samples per axis */
	moffatstruct	*moffat;	/* Array of Moffat fits to PSF */
	moffatstruct	*pfmoffat;	/* Array of pixel-free Moffat fits to PSF */
	float		moffat_fwhm_min;
	float		moffat_fwhm;	/* Central Moffat FWHM */
	float		moffat_fwhm_max;
	float		moffat_fwhm_wcs_min;
	float		moffat_fwhm_wcs;	/* Average Moffat FWHM in arcsec*/
	float		moffat_fwhm_wcs_max;
	float		moffat_ellipticity_min;
	float		moffat_ellipticity;	/* Central Moffat ellipticity */
	float		moffat_ellipticity_max;
	float		moffat_ellipticity1_min;
	float		moffat_ellipticity1;	/* Central Moffat e1 */
	float		moffat_ellipticity1_max;
	float		moffat_ellipticity2_min;
	float		moffat_ellipticity2;	/* Central Moffat e2 */
	float		moffat_ellipticity2_max;
	float		moffat_beta_min;
	float		moffat_beta;	/* Central Moffat beta */
	float		moffat_beta_max;
	float		moffat_residuals_min;
	float		moffat_residuals;/* Central Moffat residuals */
	float		moffat_residuals_max;
	float		moffat_score_min;
	float		moffat_score;	/* Central pixel-free Moffat score */
	float		moffat_score_max;
	float		pfmoffat_fwhm_min;
	float		pfmoffat_fwhm;	/* Central pixel-free Moffat FWHM */
	float		pfmoffat_fwhm_max;
	float		pfmoffat_fwhm_wcs_min;
	float		pfmoffat_fwhm_wcs; /* Average pixel-free Moffat FWHM in arcsec*/
	float		pfmoffat_fwhm_wcs_max;
	float		pfmoffat_ellipticity_min;
	float		pfmoffat_ellipticity;	/* Central pix-free Moffat ellipticity*/
	float		pfmoffat_ellipticity_max;
	float		pfmoffat_ellipticity1_min;
	float		pfmoffat_ellipticity1;	/* Central pix-free Moffat e1 */
	float		pfmoffat_ellipticity1_max;
	float		pfmoffat_ellipticity2_min;
	float		pfmoffat_ellipticity2;	/* Central pix-free Moffat e2 */
	float		pfmoffat_ellipticity2_max;
	float		pfmoffat_beta_min;
	float		pfmoffat_beta;	/* Central pixel-free Moffat beta */
	float		pfmoffat_beta_max;
	float		pfmoffat_residuals_min;
	float		pfmoffat_residuals;/* Central pixel-free Moffat residuals */
	float		pfmoffat_residuals_max;
	float		sym_residuals_min;
	float		sym_residuals;/* Symmetry residuals */
	float		sym_residuals_max;
	float		noiseqarea_min;
	float		noiseqarea;	/* Noise equivalent area */
	float		noiseqarea_max;
	float		pixscale_wcs_min;
	float		pixscale_wcs;	/* Average pixel scale in arcsec */
	float		pixscale_wcs_max;
	float		*homo_kernel;		/* PSF homogenization kernel */
	double	homopsf_params[2];	/* Idealised Moffat PSF params*/
	int		homobasis_number;	/* nb of supersampled pixels */
}	psfstruct;

typedef struct
{
	char		**command_line;			/* Command line */
	int		ncommand_line;			/* nb of params */
	char		prefs_name[MAXCHAR];		/* prefs filename */
	char		*(incat_name[MAXFILE]);		/* Filename(s) of input cats */
	int		ncat;				/* Number of input images */
	char		psf_dir[MAXCHAR];		/* PSF output dir */
	char		psf_suffix[MAXCHAR];		/* Suffix for PSF filenames */
	int		psf_size[2], npsf_size;		/* PSF size */
	enum {NEWBASIS_NONE, NEWBASIS_PCAINDEPENDENT, NEWBASIS_PCACOMMON}
	newbasis_type;			/* Type of new basis */
	int		newbasis_number;		/* Number of PCs */
	/* Point-source sample */
	double	minsn;				/* Minimum S/N for patterns */
	double	maxellip;			/* Maximum (A-B)/(A+B) */
	double	maxvar;				/* Maximum FWHM variability */
	double	fwhmrange[2];			/* Allowed FWHM range */
	int		nfwhmrange;	       		/* nb of params */
	int		flag_mask;			/* Rejection mask on SEx FLAGS*/
	int		wflag_mask;			/* Rej. mask on FLAGS_WEIGHT */
	int		imaflag_mask;			/* Rej. mask on IMAFLAGS_ISO */
	int		nmax;				/* Max. nb of samples per set*/
	double	prof_accuracy;			/* Required PSF accuracy */
	double	psf_step;			/* Oversampling (pixels) */
	double	psf_pixsize[2];			/* Eff. pixel size (pixels) */
	int		npsf_pixsize;			/* nb of params */
	int		badpix_flag;			/* Filter bad pixels? */
	int		badpix_nmax;			/* Max number of bad pixels */
	char		photflux_key[MAXCHAR];		/* Name of phot. flux key */
	char		photfluxerr_key[MAXCHAR];	/* Name of phot. flux err. key*/
	/* Vector basis */
	basistypenum	basis_type;			/* PSF vector basis set */
	int		basis_number;			/* nb of supersampled pixels */
	char		basis_name[MAXCHAR];		/* PSF vector basis filename */
	double	basis_scale;			/* Gauss-Laguerre beta param */
	/* Differential Geometry Maps */
	int		dgeo_flag;			/* Apply diff. geom. maps? */
	/* Re-centering */
	char		*(center_key[2]);		/* Names of centering keys */
	int		ncenter_key;			/* nb of params */
	int		autoselect_flag;		/* Auto. select FWHMs ? */
	int		recenter_flag;			/* Recenter PSF-candidates? */
	/* Check-images */
	checkenum	check_type[MAXCHECK];		/* check-image types */
	int		ncheck_type;			/* nb of params */
	char		*(check_name[MAXCHECK]);	/* check-image names */
	int		ncheck_name;			/* nb of params */
	int		check_cubeflag;			/* check-images as datacubes?*/
	/* PSF variability */
	enum {VAR_NONE, VAR_SEEING}	var_type;	/* PSF variability type */
	char		*(context_name[MAXCONTEXT]);	/* Names of context-keys */
	int		ncontext_name;			/* nb of params */
	int		context_group[MAXCONTEXT];	/* Context group */
	int		ncontext_group;			/* nb of params */
	int		context_nsnap;			/* nb of snapshots / context */
	int		group_deg[MAXCONTEXT];		/* Degree for each group */
	int		ngroup_deg;			/* nb of params */
	enum	{HIDDEN_MEF_INDEPENDENT, HIDDEN_MEF_COMMON}
	hidden_mef_type;		/* Mosaic handling for hiddens*/
	enum	{STABILITY_EXPOSURE, STABILITY_SEQUENCE}
	stability_type;			/* PSF stability range*/
	enum	{PSF_MEF_INDEPENDENT, PSF_MEF_COMMON}
	psf_mef_type;			/* Mosaic handling for PSF */
	/* Homogenisation kernel vector basis */
	enum	{HOMOBASIS_NONE, HOMOBASIS_GAUSSLAGUERRE}
	homobasis_type;			/* Homo. kernel basis set */
	int		homobasis_number;		/* nb of supersampled pixels */
	char		homokernel_dir[MAXCHAR];	/* Homo. kernel output dir */
	char		homokernel_suffix[MAXCHAR];	/* Homo. kernel file suffix */
	double	homobasis_scale;		/* Gauss-Laguerre beta param */
	double	homopsf_params[2];		/* Idealised Moffat PSF params*/
	int		nhomopsf_params;		/* nb of params */
	/* Output catalogs */
	char		outcat_name[MAXCHAR];		/* Output filename */
	cattypenum	outcat_type;			/* Output catalog type */
	int		outcatpipe_flag;		/* Pipe output catalogs? */
	/* Check-plots */
	cplotenum	cplot_device[MAXCHECK];		/* check-plot format */
	int		ncplot_device;			/* nb of params */
	cplotenum	cplot_type[MAXCHECK];		/* check-plot types */
	int		ncplot_type;			/* nb of params */
	char		*(cplot_name[MAXCHECK]);	/* check-plot names */
	int		ncplot_name;			/* nb of params */
	int		cplot_flag;			/* = 0 if no check-plot */
	int		cplot_res[2];			/* X,Y check-plot resolution */
	int		ncplot_res;			/* nb of params */
	int		cplot_antialiasflag;		/* Anti-aliasing on/off */
	/* Multithreading */
	int		nthreads;			/* Number of active threads */
	/* Misc */
	enum {QUIET, NORM, LOG, FULL}	verbose_type;	/* How much it displays info */
	int		xml_flag;			/* Write XML file? */
	char		xml_name[MAXCHAR];		/* XML file name */
	char		xsl_name[MAXCHAR];		/* XSL file name (or URL) */
	char		sdate_start[12];		/* PSFEx start date */
	char		stime_start[12];		/* PSFEx start time */
	char		sdate_end[12];			/* PSFEx end date */
	char		stime_end[12];			/* PSFEx end time */
	double	time_diff;			/* Execution time */
}	prefstruct;

typedef struct wcs
{
	int		naxis;			/* Number of image axes */
	int		naxisn[NAXIS];		/* FITS NAXISx parameters */
	char		ctype[NAXIS][9];	/* FITS CTYPE strings */
	char		cunit[NAXIS][32];	/* FITS CUNIT strings */
	double	crval[NAXIS];		/* FITS CRVAL parameters */
	double	cdelt[NAXIS];		/* FITS CDELT parameters */
	double	crpix[NAXIS];		/* FITS CRPIX parameters */
	double	crder[NAXIS];		/* FITS CRDER parameters */
	double	csyer[NAXIS];		/* FITS CSYER parameters */
	double	cd[NAXIS*NAXIS];	/* FITS CD matrix */
	double	*projp;			/* FITS PV/PROJP mapping parameters */
	int		nprojp;			/* number of useful projp parameters */
	double	longpole,latpole;	/* FITS LONGPOLE and LATPOLE */
	double	wcsmin[NAXIS];		/* minimum values of WCS coords */
	double	wcsmax[NAXIS];		/* maximum values of WCS coords */
	double	wcsscale[NAXIS];	/* typical pixel scale at center */
	double	wcsscalepos[NAXIS];	/* WCS coordinates of scaling point */
	double	wcsmaxradius;		/* Maximum distance to wcsscalepos */
	int		outmin[NAXIS];		/* minimum output pixel coordinate */
	int		outmax[NAXIS];		/* maximum output pixel coordinate */
	int		lat,lng;		/* longitude and latitude axes # */
	double	r0;			/* projection "radius" */
	double	lindet;			/* Determinant of the local matrix */
	int		chirality;		/* Chirality of the CD matrix */
	double	pixscale;		/* (Local) pixel scale */
	double	ap2000,dp2000;		/* J2000 coordinates of pole */
	double	ap1950,dp1950;		/* B1950 coordinates of pole */
	double	obsdate;		/* Date of observations */
	double	equinox;		/* Equinox of observations */
	double	epoch;			/* Epoch of observations (deprec.) */
	enum {RDSYS_ICRS, RDSYS_FK5, RDSYS_FK4, RDSYS_FK4_NO_E, RDSYS_GAPPT}
	radecsys;		/* FITS RADECSYS reference frame */
	celsysenum	celsys;			/* Celestial coordinate system */
	double	celsysmat[4];		/* Equ. <=> Cel. system parameters */
	int		celsysconvflag;		/* Equ. <=> Cel. conversion needed? */
	struct wcsprm	*wcsprm;		/* WCSLIB's wcsprm structure */
	struct linprm	*lin;			/* WCSLIB's linprm structure */
	struct celprm	*cel;			/* WCSLIB's celprm structure */
	struct prjprm *prj;			/* WCSLIB's prjprm structure */
	struct tnxaxis *tnx_latcor;		/* IRAF's TNX latitude corrections */
	struct tnxaxis *tnx_lngcor;		/* IRAF's TNX longitude corrections */
	struct poly	*inv_x;			/* Proj. correction polynom in x */
	struct poly	*inv_y;			/* Proj. correction polynom in y */
}	wcsstruct;

typedef struct outsample
{
	int		detindex;		/* Detection index */
	short		extindex;		/* Extension index */
	int		catindex;		/* Catalog index */
	int		badflag;		/* Rejection flags (0 = accepted) */
	double	context[MAXCONTEXT];	/* Context vector */
	int		ncontext;		/* Number of contexts */
	float		norm;			/* Normalisation flux */
	float		fwhm;			/* FWHM (derived from input catalog) */
	float		ellip;			/* Ellipticity */
	float		snr;			/* Signal-to-Noise ratio */
	double	x,y;			/* x,y position estimate in frame */
	float		dx,dy;			/* x,y shift / vignet center */
	float		chi2;			/* Chi2 of the fit */
	float		modresi;		/* Residual index */
}	outsamplestruct;

typedef struct outcat
{
	outsamplestruct	outsample;	/* Current output line */
	FILE			*ascfile;	/* Output ASCII file (if needed) */
	tabstruct		*objtab;	/* Output object table */
	keystruct		*objkeys;  	/* List of output catalog keys */
	char			*buf;		/* Line buffer */
	int			ncontext;	/* Number of contexts */
}	outcatstruct;



typedef struct field
{
	char		catname[MAXCHAR];	/* Input catalog filename */
	char		*rcatname;		/* "Reduced" catalog name */
	char		rtcatname[MAXCHAR];	/* "Reduced", no trail catalog name */
	char		ident[MAXCHAR];		/* Field identifier (read from FITS) */
	int		next;			/* Number of extensions */
	int		ndet;			/* Number of detections (info only) */
	psfstruct	**psf;			/* Array of PSFs */
	wcsstruct	**wcs;			/* Array of WCS structures */
	setstruct	*set;			/* Catalog or array of catalogs */
	catstruct	**ccat;			/* Pointers to check-image files */
	double	meanwcspos[NAXIS];	/* Mean pixel coordinate */
	double	meanwcsscale[NAXIS];	/* Mean pixel scale */
	double	maxradius;		/* Maximum radius */
	int		**lcount;		/* Count detections that are loaded */
	int		**acount;		/* Count detections that are accepted */
	int		**count;		/* Count detections in stats */
	double	**modchi2;		/* Sum of chi2's per image area */
	double	**modresi;		/* Sum of res. indices per image area */
}	fieldstruct;

//devicestruct not here cause used only in cplot

typedef struct context
{
	char		**name;			/* Context names */
	int		*group;			/* Context groups */
	int		*pcflag;		/* Flags PC contexts */
	int		ncontext;		/* Total number of contexts */
	int		*degree;		/* Group degrees */
	int		ngroup;			/* Number of context groups */
	double	*pc;			/* PC components */
	int		npc;			/* Number of PC components */
}	contextstruct;
#endif
