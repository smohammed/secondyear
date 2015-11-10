# Plot WD map

from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

wds = Table.read('GAIS_vphas_2mass_wd.txt', format='ascii')

glval, gbval, numwd = [], [], []

for i in range(0, 360):
    for j in range(-10, 10):
        if len(np.where((wds['gl_galex'] > i) & (wds['gl_galex'] < i+1) & (wds['gb_galex'] > j) & (wds['gb_galex'] < j+1))[0]) > 0.:
            glval.append(i)
            gbval.append(j)
            numwd.append(len(np.where((wds['gl_galex'] > i) & (wds['gl_galex'] < i+1) & (wds['gb_galex'] > j) & (wds['gb_galex'] < j+1))[0]))

glval = np.array(glval)
gbval = np.array(gbval)
numwd = np.array(numwd)

gl = wds['gl_galex'].tolist()
gb = wds['gb_galex'].tolist()
dresult = query(gl,gb,coordsys='gal')

dist = 10**(1.+np.array(dresult['distmod'])/5.)/1000.

for i in range(len(gl)):
    plt.plot(dist,dresult['best'][i],c=gl[i])
