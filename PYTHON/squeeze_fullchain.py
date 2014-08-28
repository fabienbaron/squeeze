import numpy as np
import math
from astropy.io import fits
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

fullchain_filename = "../pipo_fullchain.fits.gz"
hdu_list = fits.open(fullchain_filename)
header = hdu_list[0].header
image_data = hdu_list[0].data
sz = image_data.shape

nthreads = header['nthreads']
niter = header['niter']
nchanr   = header['nchanr']
nelements = header['elements']
axis_len = sz[3] 
snplots = math.ceil(np.sqrt(niter))
print('Nthreads: ', nthreads)
print('Niter: ', niter)
print('Nchanr: ', nchanr)
print('Nelements: ', nelements)
print('Image width: ', axis_len)
print('Image grid: ', snplots)
plt.set_cmap('hot') 

fig = plt.figure(1, (20., 20.), frameon = False)

#ax = plt.subplot(111)
#plt.axis('off')
#plt.setp(ax, 'frame_on', False)
#ax.set_xticks([])
#ax.set_yticks([])
#ax.grid('off')


grid = ImageGrid(fig, 111, nrows_ncols = (snplots, snplots), axes_pad=0.0, cbar_mode=None, share_all=True)

fig = plt.figure(1, (20., 20.), frameon = False)
ax = fig.add_axes([0, 0, axis_len, axis_len])
ax.axis('off')

grid = ImageGrid(fig, 111, nrows_ncols = (snplots, snplots), axes_pad=0.0)

for i in range(snplots*snplots):
    if i < niter:
        grid[i].imshow(np.sqrt(image_data[0,i,0,:,:].reshape(axis_len, axis_len)))
    else: 
        grid[i].imshow(np.zeros(shape=(axis_len, axis_len)))


plt.show()
