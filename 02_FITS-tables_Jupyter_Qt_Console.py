# Jupyter Qt Console:
# https://qtconsole.readthedocs.io/en/stable/installation.html

# Explicação "%matplotlib inline" e exemplo:
# http://nbviewer.jupyter.org/github/ipython/ipython/blob/1.x/examples/notebooks/Part%203%20-%20Plotting%20with%20Matplotlib.ipynb

# Astropy Exemplo:
# http://www.astropy.org/astropy-tutorials/FITS-tables.html

import numpy as np

from astropy.io import fits

import matplotlib.pyplot as plt

%matplotlib inline


# The following line is needed to download the example FITS files used in this tutorial.

from astropy.utils.data import download_file

#

from astropy.io import fits


# link para baixar o chandra_events.fits: https://astropy.stsci.edu/data/tutorials/FITS-tables/

# Tamanho do arquivo: 30.3Mb

# event_filename = ('C:\chandra_events.fits')


event_filename = download_file( 'http://data.astropy.org/tutorials/FITS-tables/chandra_events.fits', cache=True )

# Pasta do cache no Windows: Filename: C:\Users\Teste\.astropy\cache\download\py3\

hdu_list = fits.open(event_filename, memmap=True)

hdu_list.info()

print(hdu_list[1].columns)

from astropy.table import Table

evt_data = Table(hdu_list[1].data)

evt_data

NBINS = 500
energy_hist = plt.hist(evt_data['energy'], NBINS)

ii = np.in1d(evt_data['ccd_id'], [0, 1, 2, 3])
np.sum(ii)

NBINS = (100,100)

img_zero, yedges, xedges = np.histogram2d(evt_data['x'][ii], evt_data['y'][ii], NBINS)

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

plt.imshow(img_zero, extent=extent, interpolation='nearest', cmap='gist_yarg', origin='lower')

plt.xlabel('x')
plt.ylabel('y')

# To see more color maps
# http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps

from matplotlib.colors import LogNorm

NBINS = (100,100)
img_zero_mpl = plt.hist2d(evt_data['x'][ii], evt_data['y'][ii], NBINS, cmap='viridis', norm=LogNorm())

cbar = plt.colorbar(ticks=[1.0,3.0,6.0])
cbar.ax.set_yticklabels(['1','3','6'])

plt.xlabel('x')
plt.ylabel('y')
