"""
Read the psfex generated representation of the psf
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
hdu_list = fits.open('test.psf')
print(hdu_list.info())

data = hdu_list[1].data

hdu_list.close()
print(data)


nrows = 3
ncols = 2
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
ax = ax.ravel()
for i in range(nrows*ncols):
    ax[i].imshow(data[0][0][i], origin='lower', cmap='viridis')
plt.show()

fdata = np.zeros((25,25))
for i in range(6):
    fdata += data[0][0][i]
plt.imshow(fdata, origin='lower', cmap='viridis')
plt.colorbar()
plt.show()
plt.imshow(data[0][0][0], origin='lower', cmap='viridis')
plt.show()
