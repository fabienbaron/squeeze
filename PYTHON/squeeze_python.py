import numpy as np

import matplotlib 
import matplotlib.pyplot as plt
from astropy.io import fits


thread_filename = "../pipo.fits"
hdu_list = fits.open(thread_filename)
image_data = hdu_list[0].data
sz = image_data.shape
image_data = image_data.reshape(sz[1], sz[2])
plt.ioff()
plt.imshow(image_data, cmap='hot', interpolation='none')
plt.colorbar()
plt.show()

#snplots = ceil(sqrt(niter))
#gs = gridspec.GridSpec(snplots, snplots, wspace=0.0)
