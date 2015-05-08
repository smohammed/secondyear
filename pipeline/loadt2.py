import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky

# Load tables, make cuts
tycho = fits.open('../tychobstarmatch.fits')[1].data
tcut = np.where((tycho.gb > -10) & (tycho.gb < 10) & (tycho.nuv_cps > 10.) & (tycho.nuv_cps < 100.))
tycho = tycho[tcut]
ext = fits.open('../dust/extinction.fits')[1].data
t2 = tycho[:len(ext)]

# Match stars with extinction catalog
t2gal = SkyCoord(t2.ra*u.degree,t2.dec*u.degree,frame='icrs')
extgal = SkyCoord(ext.ra*u.degree,ext.dec*u.degree,frame='icrs')
t2ind, extind, angsep, dist3d = search_around_sky(t2gal, extgal, 0.1*u.arcsec)
t2 = t2[t2ind]
ext = ext[extind]


# Bands
#nuv = t2.nuvmag
#fuv = t2.fuvmag
#BJ = t2.BJmag
#VJ = t2.VJmag

# Extinction A_lambda, from Wyder et al 2007
rnuv = 8.2
rfuv = 8.24
anuv = rnuv * ext.E_B_V_SandF
afuv = rfuv * ext.E_B_V_SandF

aV = ext.AV_SandF
aB = aV + ext.E_B_V_SandF

#nuvext = t2.nuvmag - anuv
#fuvext = t2.fuvmag - afuv

