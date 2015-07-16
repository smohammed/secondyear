import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky

# Load tables, make cuts
tycho = fits.open('../tychobstarmatch.fits')[1].data
tcut = np.where((tycho.gb > -10) & (tycho.gb < 10) & (tycho.nuv_cps > 5.) & (tycho.nuv_cps < 100.))
tycho = tycho[tcut]
ext = fits.open('../dust/extinction.fits')[1].data
t2 = tycho
photons = fits.open('../photometry12000.fits')[1].data
sdss = fits.open('../sdss_0-300.fits')[1].data
sex = fits.open('../combmaps12000/sex_total.fits')[1].data

photons = sex

# Match stars with photometry
t2gal = SkyCoord(t2.gl*u.degree, t2.gb*u.degree, frame='galactic')
photgal = SkyCoord(photons.gl*u.degree, photons.gb*u.degree, frame='galactic')
t2ind, photind, angsep, dist3d = search_around_sky(t2gal, photgal, 3.*u.arcsec)
t2 = t2[t2ind]
photons = photons[photind]

# Match stars with extinction catalog
t2icrs = SkyCoord(t2.ra*u.degree, t2.dec*u.degree, frame='icrs')
exticrs = SkyCoord(ext.ra*u.degree, ext.dec*u.degree, frame='icrs')
t2newind, extind, angsep, dist3d = search_around_sky(t2icrs, exticrs, 0.1*u.arcsec)
t2 = t2[t2newind
]ext = ext[extind]
photons = photons[t2newind]

cuts = np.unique(t2.galex_id, return_index='true')[1]
t2 = t2[cuts]
ext = ext[cuts]
photons = photons[cuts]

# Extinction A_lambda, from Wyder et al 2007
rnuv = 8.2
rfuv = 8.24
anuv = rnuv * ext.E_B_V_SandF
extcut = np.where(anuv < 10)

t2 = t2[extcut]
ext = ext[extcut]
photons = photons[extcut]

anuv = rnuv * ext.E_B_V_SandF
afuv = rfuv * ext.E_B_V_SandF
aV = ext.AV_SandF
aB = aV + ext.E_B_V_SandF

#snuv = photons.snuv * 10.

# Now match sdss catalog
t2icrs = SkyCoord(t2.ra*u.degree, t2.dec*u.degree, frame='icrs')
sdssicrs = SkyCoord(sdss.ra*u.degree, sdss.dec*u.degree, frame='icrs')
t2sdssind, sdssind, angsep, dist3d = search_around_sky(t2icrs, sdssicrs, 5.*u.arcsec)
t3 = t2[t2sdssind]
ext3 = ext[t2sdssind]
photons3 = photons[t2sdssind]
sdss3 = sdss[sdssind]

#snuv3 = photons3.snuv * 10.
