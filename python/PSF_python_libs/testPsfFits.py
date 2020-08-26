"""
generate ePSF from an input FITS file
"""

from astropy.visualization import simple_norm
from photutils import datasets
from astropy.io import fits
from astropy.utils.data import download_file

hdu_list = fits.open('/mnt/NewHDD/testPSFEx/SourceExtractor/in/x0bq0102t_c0h.fit')
print(hdu_list.info())

data = hdu_list[0].data

hdu_list.close()

import matplotlib.pyplot as plt
plt.imshow(data, cmap='gray')
plt.colorbar()
plt.show()

from photutils import find_peaks
peaks_tbl = find_peaks(data, threshold=20.)

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

from photutils import EPSFBuilder
epsf_builder = EPSFBuilder(oversampling=4, maxiters=3,
                           progress_bar=False)
epsf, fitted_stars = epsf_builder(stars)


norm = simple_norm(epsf.data, 'log', percent=99.)
plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
plt.colorbar() 
plt.show() 
