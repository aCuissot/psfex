"""
using tutorial https://photutils.readthedocs.io/en/stable/psf.html
"""

from photutils.psf import IterativelySubtractedPSFPhotometry
import numpy as np
from astropy.table import Table
from photutils.datasets import make_random_gaussians_table, make_noise_image, make_gaussian_sources_image
from matplotlib import rcParams
import matplotlib.pyplot as plt
from photutils.psf.photometry import BasicPSFPhotometry

sigma_psf = 2.0
sources = Table()
sources['flux'] = [700, 800, 700, 800]
sources['x_mean'] = [12, 17, 12, 17]
sources['y_mean'] = [15, 15, 20, 20]
sources['x_stddev'] = sigma_psf * np.ones(4)
sources['y_stddev'] = sources['x_stddev']
sources['theta'] = [0, 0, 0, 0]
sources['id'] = [1, 2, 3, 4]
tshape = (32, 32)
image = (make_gaussian_sources_image(tshape, sources) + 
         make_noise_image(tshape, distribution='poisson', mean=6., random_state=1) +
         make_noise_image(tshape, distribution='gaussian', mean=0., stddev=2, random_state=1))

plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()

from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm

bkgrms = MADStdBackgroundRMS()
std = bkgrms(image)
iraffind = IRAFStarFinder(threshold=3.5*std, fwhm=sigma_psf*gaussian_sigma_to_fwhm, minsep_fwhm=0.01,
                          roundhi=5.0 ,roundlo=-5.0, sharplo=0.0, sharphi=2.0)
daogroup = DAOGroup(2.0*sigma_psf*gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=sigma_psf)

photometry = IterativelySubtractedPSFPhotometry(finder=iraffind, group_maker=daogroup, bkg_estimator=mmm_bkg,
                                                psf_model=psf_model, fitter=LevMarLSQFitter(), niters=1, fitshape=(11, 11))
result_tab = photometry(image=image)
residual_image = photometry.get_residual_image()

plt.subplot(1, 2, 1)
plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.title('Simulated data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(residual_image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.title('Residual image')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()

psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
pos = Table(names=['x_0', 'y_0'], data=[sources['x_mean'], sources['y_mean']])

photometry = BasicPSFPhotometry(group_maker=daogroup, bkg_estimator=mmm_bkg, psf_model=psf_model,
                                fitter=LevMarLSQFitter(), fitshape=(11,11))

result_tab = photometry(image=image, init_guesses=pos)
residual_image = photometry.get_residual_image()

plt.subplot(1, 2, 1)
plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.title('Simulated data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(residual_image, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.title('Residual image')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()

gaussian_prf = IntegratedGaussianPRF()
gaussian_prf.sigma.fixed = False
gaussian_prf.sigma.value = 2.05

sources = Table()
sources['flux'] = [10000, 1000]
sources['x_mean'] = [18, 9]
sources['y_mean'] = [17, 21]
sources['x_stddev'] = [2] * 2
sources['y_stddev'] = sources['x_stddev']
sources['theta'] = [0] * 2
tshape = (32, 32)

from matplotlib.colors import LogNorm

image = (make_gaussian_sources_image(tshape, sources) + 
         make_noise_image(tshape, distribution='poisson', mean=6., random_state=1) +
         make_noise_image(tshape, distribution='gaussian', mean=0., stddev=2, random_state=1))

vmin, vmax = np.percentile(image, [5, 95])

plt.imshow(image, cmap='viridis', aspect=1, interpolation='nearest', 
           origin='lower', norm=LogNorm(vmin=vmin, vmax=vmax))
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()

daogroup = DAOGroup(crit_separation=8)
mmm_bkg = MMMBackground()
iraffind = IRAFStarFinder(threshold=2.5*mmm_bkg(image), fwhm=4.5)
fitter = LevMarLSQFitter()
gaussian_prf = IntegratedGaussianPRF(sigma=2.05)
gaussian_prf.sigma.fixed = False
itr_phot_obj = IterativelySubtractedPSFPhotometry(finder=iraffind, group_maker=daogroup, bkg_estimator=mmm_bkg,
                                                  psf_model=gaussian_prf, fitter=fitter, fitshape=(11, 11), niters=2)

phot_results = itr_phot_obj(image)

print(phot_results['id', 'group_id', 'iter_detected', 'x_0', 'y_0', 'flux_0'])
print(phot_results['sigma_0', 'sigma_fit', 'x_fit', 'y_fit', 'flux_fit'])

plt.imshow(itr_phot_obj.get_residual_image(), cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()
