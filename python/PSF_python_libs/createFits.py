"""
Generate a sky image with some stars and save it at the FITS format
"""

from astropy.visualization import simple_norm
from photutils import datasets
import matplotlib.pyplot as plt

hdu = datasets.load_simulated_hst_star_image()
data = hdu.data
print(data.shape)
from photutils.datasets import make_noise_image
data +=  make_noise_image(data.shape, distribution='gaussian', mean=10.,
                          stddev=5., random_state=12345)
norm = simple_norm(data, 'sqrt', percent=99.)
plt.imshow(data, norm=norm, origin='lower', cmap='gray')
plt.show()

from astropy.io import fits
new_hdul = fits.HDUList()
img_hdu = fits.ImageHDU(data=data)
img_hdu.update_header()
new_hdul.append(img_hdu)
new_hdul.update_extend()

new_hdul.writeto('/tmp/test.fit')
