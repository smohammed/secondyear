from intphotometry import phot
#from loadt2 import *
import itertools
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky

'''
# Get photometry
sgl = []
sgb = []
table = []

for i in range(0, 360, 10):
	gl, gb, t = phot(i, i+10, 0, 10)
	sgl.append(gl)
	sgb.append(gb)
	table.append(t)
for i in range(0, 360, 10):
    gl, gb, t = phot(i, i+10, -10, 0)
    sgl.append(gl)
    sgb.append(gb)
    table.append(t)
sglt = np.array(list(itertools.chain(*sgl)))
sgbt = np.array(list(itertools.chain(*sgb)))
snuv = np.array(list(itertools.chain(*table)))

from numpy import inf
snuv[snuv == inf] = 0
q = np.isnan(snuv)
snuv[q] = 0
cuts = np.where(snuv > 0)
sglt = sglt[cuts]
sgbt = sgbt[cuts]
snuv = -2.5*np.log10(snuv[cuts]) + 20.08
'''


phot = fits.open('../phottable.fits')[1].data

t2sgal = SkyCoord(t2.gl*u.degree,t2.gb*u.degree,frame='galactic')
sgal = SkyCoord(phot.gl*u.degree,phot.gb*u.degree,frame='galactic')

sind, t2sind, angseps, dist3ds = search_around_sky(sgal, t2sgal, 0.1*u.arcsec)

t2s = t2[t2sind]
phot = phot[sind]
pgl = phot.gl
pgb = phot.gb
nuv = phot.nuv

