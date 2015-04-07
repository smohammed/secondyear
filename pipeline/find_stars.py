import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky
matplotlib.rcParams['figure.figsize'] = 16,8
matplotlib.rcParams['font.size'] = 20

###################################################################
# Load star values
###################################################################
'''
x = np.loadtxt('x_stars.txt',dtype='string').astype(float)
y = np.loadtxt('y_stars.txt',dtype='string').astype(float)
x = x * 10./12000.
y = y * 10./12000. - 10. 
'''
star_cand = fits.open('find_stars.fits')[1].data
sgl = star_cand.gl.tolist()
sgb = star_cand.gb.tolist()
s_cgal = SkyCoord(sgl*u.degree, sgb*u.degree, frame='galactic')

###################################################################
# Load bstars
###################################################################
bstar = fits.open('bstar.fits')[1].data
bcut1 = np.where((bstar.nuv_cps > 10.) & (bstar.nuv_cps < 100.))
bstar = bstar[bcut1]
bgal = SkyCoord(bstar.ra*u.degree, bstar.dec*u.degree, frame='icrs').galactic
bcut2 = np.where((bgal.l.degree > 0) & (bgal.l.degree < 10) & (bgal.b.degree > -10.) & (bgal.b.degree < 10.))
bstar2 = bstar[bcut2]
bgal2 = SkyCoord(bstar2.ra*u.degree, bstar2.dec*u.degree, frame='icrs').galactic


###################################################################
# Match stars with detections
###################################################################
bstarind, star_candind, starsep, stardist = search_around_sky(bgal2,s_cgal,1.*u.arcmin)





bstar2 = bstar2[np.unique(bstarind)]
bstarname = bstar2.galex_id
bgal2 = SkyCoord(bstar2.ra*u.degree, bstar2.dec*u.degree, frame='icrs').galactic
bgl = bgal2.l.degree
bgb = bgal2.b.degree


bstarnumber = []
for i in range(len(bstarind)-1):
	if bstarind[i] != bstarind[i+1]:
		bstarnumber.append(i) 		# Do this so we can match indices with photind

