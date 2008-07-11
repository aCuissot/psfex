 /*
 				makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Main program.
*
*	Last modify:	08/07/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"check.h"
#include	"context.h"
#include	"cplot.h"
#include	"diagnostic.h"
#include	"field.h"
#include	"homo.h"
#include	"pca.h"
#include	"prefs.h"
#include	"psf.h"
#include	"sample.h"
#include	"xml.h"

psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context);
time_t		thetime, thetime2;

/********************************** makeit ***********************************/
/*
*/
void	makeit(void)

  {
   fieldstruct		**fields;
   psfstruct		**cpsf,
			*psf;
   setstruct		*set, *set2;
   contextstruct	*context, *fullcontext;
   struct tm		*tm;
   static char		str[MAXCHAR];
   char			**incatnames,
			*pstr;
   float		*psfsteps, *basis,
			psfstep;
   int			c,i, ncat, ext, next, nmed, nbasis;

/* Install error logging */
  error_installfunc(write_error);

  incatnames = prefs.incat_name;
  ncat = prefs.ncat;

/* Processing start date and time */
  thetime = time(NULL);
  tm = localtime(&thetime);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT,
	"----- %s %s started on %s at %s with %d thread%s\n\n",
		BANNER,
		MYVERSION,
		prefs.sdate_start,
		prefs.stime_start,
		prefs.nthreads,
		prefs.nthreads>1? "s":"");


/* End here if no filename has been provided */
  if (!ncat)
    {
/*-- Processing end date and time */
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);

/*-- Write XML */
    if (prefs.xml_flag)
      {
      init_xml(0);
      write_xml(prefs.xml_name);
      end_xml();
      }
    return;
    }

/* Create an array of PSFs (one PSF for each extension) */
  QMALLOC(fields, fieldstruct *, ncat);
  for (c=0; c<ncat; c++)
    fields[c] = field_init(incatnames[c]);

  next = fields[0]->next;

  if (prefs.xml_flag)
    init_xml(ncat);  

  psfstep = prefs.psf_step;
  psfsteps = NULL;
  nbasis = 0;
  basis = NULL;

/* Initialize context */
  context = context_init(prefs.context_name, prefs.context_group,
		prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
		CONTEXT_REMOVEPC);
  fullcontext = context->npc?
		  context_init(prefs.context_name, prefs.context_group,
		    prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
		    CONTEXT_KEEPPC)
		: context;

  if (prefs.newbasis_type==NEWBASIS_PCAMULTI && (!psfstep))
    {
/* A first run through all samples only to derive a common pixel step */
    QMALLOC(psfsteps, float, next);
    for (ext=0 ; ext<next; ext++)
      {
      set = load_samples(incatnames, ncat, ext, next, context);
      psfsteps[ext] = (float)(psfstep? psfstep : (set->fwhm/2.35)*0.5);
      end_set(set);
      }
    }
  else if (prefs.newbasis_type==NEWBASIS_PCASINGLE)
    {
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT,
	"----- First pass for Principal Component Analysis (single): "
	"%dx%d PSFs required\n", ncat, next);
    if (!prefs.psf_step)
      {
      QMALLOC(psfsteps, float, next);
      for (ext=0 ; ext<next; ext++)
        {
        set = load_samples(incatnames, ncat, ext, next, context);
        psfsteps[ext] = (float)((set->fwhm/2.35)*0.5);
        end_set(set);
        }
      psfstep = fast_median(psfsteps, next);
      }
/*-- Derive a new common PCA basis for all extensions */
    QMALLOC(cpsf, psfstruct *, ncat*next);
    for (ext=0 ; ext<next; ext++)
      for (c=0; c<ncat; c++)
        {
        set = load_samples(&incatnames[c], 1, ext, next, context);
        sprintf(str, "Computing PSF model for catalog %s...", incatnames[c]);
        NFPRINTF(OUTPUT, str);
        cpsf[c+ext*ncat] = make_psf(set, psfstep, NULL, 0, context);
        end_set(set);
        }
    nbasis = prefs.newbasis_number;
    basis = pca_onsnaps(cpsf, ncat*next, nbasis);
    for (i=0 ; i<ncat*next; i++)
      psf_end(cpsf[i]);
    free(cpsf);
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT,
	"----- Second pass for Principal Component Analysis (single):\n\n");
    }

  QIPRINTF(OUTPUT,
        " extension accepted/total sampling chi2/dof FWHM(pix) ellip."
	" residuals asymmetry");
  for (ext=0 ; ext<next; ext++)
    {
    if (prefs.newbasis_type == NEWBASIS_PCAMULTI)
/*---- Derive a new PCA basis for each extension */
      {
      psfstep = prefs.psf_step? prefs.psf_step : psfsteps[ext];
      QMALLOC(cpsf, psfstruct *, ncat);
      for (c=0; c<ncat; c++)
        {
        set = load_samples(&incatnames[c], 1, ext, next, context);
        sprintf(str, "Computing PSF model for catalog %s...", incatnames[c]);
        NFPRINTF(OUTPUT, str);
        cpsf[c] = make_psf(set, psfstep, NULL, 0, context);
        end_set(set);
        }
      nbasis = prefs.newbasis_number;
      basis = pca_onsnaps(cpsf, ncat, nbasis);
      for (c=0 ; c<ncat; c++)
        psf_end(cpsf[c]);
      free(cpsf);
      }

    if (context->npc)
/*---- Derive principal components of PSF components */
      {
      QMALLOC(cpsf, psfstruct *, ncat);
      for (c=0; c<ncat; c++)
        {
        set = load_samples(&incatnames[c], 1, ext, next, context);
        if (prefs.newbasis_type == NEWBASIS_NONE && !psfstep)
          psfstep = (float)((set->fwhm/2.35)*0.5);
        sprintf(str, "Computing PSF model for catalog %s...", incatnames[c]);
        NFPRINTF(OUTPUT, str);
        cpsf[c] = make_psf(set, psfstep, basis, nbasis, context);
        end_set(set);
        }
      free(fullcontext->pc);
      fullcontext->pc = pca_oncomps(cpsf, ncat, context->npc);
      for (c=0 ; c<ncat; c++)
        psf_end(cpsf[c]);
      free(cpsf);
      }

/*-- Load all the samples */
    set = load_samples(incatnames, ncat, ext, next, fullcontext);
    if (prefs.newbasis_type == NEWBASIS_NONE && !psfstep)
      psfstep = (float)((set->fwhm/2.35)*0.5);
    psf = make_psf(set, psfstep, basis, nbasis, fullcontext);
    end_set(set);
    NFPRINTF(OUTPUT, "Computing final model...");
    context_apply(fullcontext, psf, fields, ext, ncat);
    nmed = psf->nmed;
    psf_end(psf);
    for (c=0; c<ncat; c++)
      {
      psf = fields[c]->psf[ext];
/*---- Check PSF with individual datasets */
      set2 = load_samples(&incatnames[c], 1, ext, next, context);
      psf->samples_loaded = set2->nsample;
      if (set2->nsample>1)
        {
/*------ Remove bad PSF candidates */
        psf_clean(psf, set2, prefs.prof_accuracy);
        psf->chi2 = set2->nsample? psf_chi2(psf, set2) : 0.0;
        }
      psf->samples_accepted = set2->nsample;
      fields[c]->set = set2;

/*---- Compute diagnostics */
      NFPRINTF(OUTPUT,"Computing diagnostics...");
      psf_diagnostic(psf);
      QPRINTF(OUTPUT, "[%3d/%-3d]     %5d/%-5d  %6.2f    %6.2f %6.2f     %4.2f"
	"    %5.2f     %5.2f\n",
	ext+1, next,
	psf->samples_accepted, psf->samples_loaded,
	psf->pixstep,
	psf->chi2,
	sqrt(psf->moffat[nmed].fwhm_min*psf->moffat[nmed].fwhm_max),
	(psf->moffat[nmed].fwhm_max-psf->moffat[nmed].fwhm_min)
	/ (psf->moffat[nmed].fwhm_max+psf->moffat[nmed].fwhm_min),
	psf->moffat[nmed].residuals, psf->moffat[nmed].symresiduals);
/*---- Save "Check-images" */
      for (i=0; i<prefs.ncheck_type; i++)
        if (prefs.check_type[i])
          {
          sprintf(str, "Saving CHECK-image #%d...", i+1);
          NFPRINTF(OUTPUT, str);
          check_write(fields[c], prefs.check_name[i], prefs.check_type[i],
		ext, next, prefs.check_cubeflag);
          }

/*---- Free memory */
      end_set(set2);
      }
    }

  free(psfsteps);

/* Save result */
  NFPRINTF(OUTPUT,"Saving the PSF descriptions...");
  for (c=0; c<ncat; c++)
    {
/*-- Create a file name with a "PSF" extension */
    strcpy(str, incatnames[c]);
    if (!(pstr = strrchr(str, '.')))
      pstr = str+strlen(str);
    sprintf(pstr, "%s", prefs.psf_suffix);
    field_psfsave(fields[c], str);
/* Create homogenisation kernels */
    if (prefs.homobasis_type != HOMOBASIS_NONE)
      {
      NFPRINTF(OUTPUT, "Computing homogenisation kernel...");
      strcpy(str, incatnames[c]);
      if (!(pstr = strrchr(str, '.')))
        pstr = str+strlen(str);
        sprintf(pstr, "%s", prefs.homokernel_suffix);
      for (ext=0; ext<next; ext++)
        psf_homo(fields[c]->psf[ext], str, prefs.homopsf_params,
		prefs.homobasis_number, prefs.homobasis_scale, ext, next);
      }
#ifdef HAVE_PLPLOT
/* Plot diagnostic maps for all catalogs */
    cplot_ellipticity (fields[c]);
    cplot_fwhm(fields[c]);
#endif
/*-- Update XML */
    if (prefs.xml_flag)
      update_xml(fields[c]);
    }

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = difftime(thetime2, thetime);

/* Write XML */
  if (prefs.xml_flag)
    {
    NFPRINTF(OUTPUT, "Writing XML file...");
    write_xml(prefs.xml_name);
    end_xml();
    }

/* Free memory */
  for (c=0; c<ncat; c++)
    field_end(fields[c]);

  if (context->npc)
    context_end(fullcontext);   
  context_end(context);   

  return;
  }


/****** make_psf *************************************************************
PROTO	psfstruct *make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context)
PURPOSE	Make PSFs from a set of FITS binary catalogs.
INPUT	Pointer to a sample set,
	PSF sampling step,
	Pointer to basis image vectors,
	Number of basis vectors,
	Pointer to context structure.
OUTPUT  Pointer to the PSF structure.
NOTES   Diagnostics are computed only if diagflag != 0.
AUTHOR  E. Bertin (IAP)
VERSION 26/03/2008
 ***/
psfstruct	*make_psf(setstruct *set, float psfstep,
			float *basis, int nbasis, contextstruct *context)
  {
   psfstruct		*psf;
   basistypenum		basistype;

  NFPRINTF(OUTPUT,"Initializing PSF modules...");
  psf = psf_init(context, prefs.psf_size, psfstep, set->nsample);

  psf->samples_loaded = set->nsample;
  psf->fwhm = set->fwhm;
  
/* Make the basic PSF-model (1st pass) */
  NFPRINTF(OUTPUT,"Modeling the PSF.");
  psf_make(psf, set, 0.2);
  if (basis && nbasis)
    {
    QMEMCPY(basis, psf->basis, float, nbasis*psf->size[0]*psf->size[1]);
    psf->nbasis = nbasis;
    }
  else
    {
    NFPRINTF(OUTPUT,"Generating the PSF model...");
    basistype = prefs.basis_type;
    if (basistype==BASIS_PIXEL_AUTO)
      basistype = (psf->fwhm < PSF_AUTO_FWHM)? BASIS_PIXEL : BASIS_NONE;
    psf_makebasis(psf, set, basistype, prefs.basis_number);
    }
  psf_refine(psf, set);

/* Remove bad PSF candidates */
  if (set->nsample>1)
    {
    psf_clean(psf, set, 0.2);

/*-- Make the basic PSF-model (2nd pass) */
    NFPRINTF(OUTPUT,"Modeling the PSF...");
    psf_make(psf, set, 0.1);
    psf_refine(psf, set);
    }

/* Remove bad PSF candidates */
  if (set->nsample>1)
    {
    psf_clean(psf, set, 0.1);

/*-- Make the basic PSF-model (3rd pass) */
    NFPRINTF(OUTPUT,"Modeling the PSF...");
    psf_make(psf, set, 0.05);
    psf_refine(psf, set);
    }

/* Remove bad PSF candidates */
  if (set->nsample>1)
    psf_clean(psf, set, 0.05);

  psf->samples_accepted = set->nsample;

/* Refine the PSF-model */
  psf_make(psf, set, prefs.prof_accuracy);
  psf_refine(psf, set);

/* Clip the PSF-model */
  psf_clip(psf);

/*-- Just check the Chi2 */
  psf->chi2 = set->nsample? psf_chi2(psf, set) : 0.0;

  return psf;
  }


/****** write_error ********************************************************
PROTO	void    write_error(char *msg1, char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/02/2007
 ***/
void	write_error(char *msg1, char *msg2)
  {
   char	error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);
  end_xml();

  return;
  }


