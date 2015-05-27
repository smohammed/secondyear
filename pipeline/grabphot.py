from intphotometry import phot
import itertools
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
#import matplotlib.cm as cm
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# Get photometry
sgl = []
sgb = []
table = []
'''
for i in range(0, 140, 20):
    print 'gl = ' + str(i)
    gl, gb, t = phot(i, i+20, -10, 10)
    sgl.append(gl)
    sgb.append(gb)
    table.append(t)
for i in range(160, 240, 20):
    print 'gl = ' + str(i)
    gl, gb, t = phot(i, i+20, -10, 10)
    sgl.append(gl)
    sgb.append(gb)
    table.append(t)
'''
for i in range(280, 360, 20):
    print 'gl = ' + str(i)
    gl, gb, t = phot(i, i+20, -10, 10)
    sgl.append(gl)
    sgb.append(gb)
    table.append(t)

sglt = np.array(list(itertools.chain(*sgl)))
sgbt = np.array(list(itertools.chain(*sgb)))
snuv = np.array(list(itertools.chain(*table)))

from numpy import inf
snuv[snuv == inf] = 0
snuv[snuv == -inf] = 0
snuv[np.isnan(snuv)] = 0
cuts = np.where(snuv > 0)
sglt = sglt[cuts]
sgbt = sgbt[cuts]
snuv = -2.5*np.log10(snuv[cuts]) + 20.08

col1 = fits.Column(name='gl', format='E', array=sglt)
col2 = fits.Column(name='gb', format='E', array=sgbt)
col3 = fits.Column(name='snuv', format='E', array=snuv)
new_cols = fits.ColDefs([col1, col2, col3])
hdu = fits.BinTableHDU.from_columns(new_cols)
hdu.writeto('../photometry_12000_260-340.fits')

plt.scatter(sglt, sgbt)
plt.show()


'''
photons = fits.open('../phottable.fits')[1].data

t2sgal = SkyCoord(t2.gl*u.degree,t2.gb*u.degree,frame='galactic')
sgal = SkyCoord(photons.gl*u.degree,photons.gb*u.degree,frame='galactic')

sind, t2sind, angseps, dist3ds = search_around_sky(sgal, t2sgal, 0.1*u.arcsec)

t2s = t2[t2sind] # Tycho 2 stars
ext = ext[t2sind] # Extinction
photons = photons[sind] # Photometry
pgl = photons.gl
pgb = photons.gb
nuv = photons.nuv

BJ = t2s.BJmag
VJ = t2s.VJmag
'''