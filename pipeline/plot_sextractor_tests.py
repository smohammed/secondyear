from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import os
import matplotlib
import matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20


img = fits.open('../combmaps12000/0_im1_rand_deadtime.fits')[0].data

os.system('sex ../combmaps12000/0_im1_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_0_im1_rand_deadtime_test.fits -CHECKIMAGE_NAME ../combmaps12000/background_0_im1_rand_deadtime_test.fits')

im1xmin, im1xmax, im1ymin, im1ymax = 0, 4400, 0, 11999

im1sex = fits.open('../combmaps12000/sex_0_im1_rand_deadtime_test.fits')[1].data
data = im1sex
xfac = im1xmin
yfac = im1ymin
dgl = (data.X_IMAGE+xfac) * 20./12000. + 0.
dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08

#new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
#hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
#hdu.writeto('../combmaps12000/sex_0_im1_rand_deadtime_edit.fits')


plt.imshow(img, vmin=0, vmax=0.1, origin='lower', interpolation='nearest', aspect='auto', cmap=cm.gray)
plt.scatter(data.X_IMAGE, data.Y_IMAGE, facecolor='none', edgecolor='red', s=20)
plt.show()
