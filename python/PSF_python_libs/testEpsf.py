"""
test the photutils epsf
"""

from astropy.visualization import simple_norm
from photutils import datasets
import matplotlib.pyplot as plt
import numpy as np


hdu = datasets.load_simulated_hst_star_image()
data = hdu.data
print(data.shape)
from photutils.datasets import make_noise_image
data +=  make_noise_image(data.shape, distribution='gaussian', mean=10.,
                          stddev=5., random_state=12345)
norm = simple_norm(data, 'sqrt', percent=99.)
plt.imshow(data, norm=norm, origin='lower', cmap='gray')
plt.show()

from photutils import find_peaks
peaks_tbl = find_peaks(data, threshold=500.)
peaks_tbl['peak_value'].info.format = '%.8g'  # for consistent table output
print(peaks_tbl)
  
size = 25
hsize = (size - 1) / 2
x = peaks_tbl['x_peak']
y = peaks_tbl['y_peak']
mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) &
        (y > hsize) & (y < (data.shape[0] -1 - hsize)))

from astropy.table import Table
stars_tbl = Table()
stars_tbl['x'] = x[mask]
stars_tbl['y'] = y[mask]

from astropy.stats import sigma_clipped_stats
mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)
data -= median_val

from astropy.nddata import NDData
nddata = NDData(data=data)

from photutils.psf import extract_stars
stars = extract_stars(nddata, stars_tbl, size=25)

import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
nrows = 5
ncols = 5
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()
for i in range(nrows*ncols):
    norm = simple_norm(stars[i], 'log', percent=99.)
    ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
plt.show()

from photutils import EPSFBuilder
epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=False)
epsf, fitted_stars = epsf_builder(stars)

norm = simple_norm(epsf.data, 'log', percent=99.)
plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')

plt.colorbar() 
plt.show()
