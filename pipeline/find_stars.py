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
x = np.loadtxt('x_stars_flux5.txt',dtype='string').astype(float)
y = np.loadtxt('y_stars_flux5.txt',dtype='string').astype(float)
flux = np.loadtxt('flux_stars_flux5.txt',dtype='string').astype(float)
x = x * 10./12000.
y = y * 10./12000. - 10.
hdu = fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='gl', format='E', array=x),fits.Column(name='gb', format='E', array=y),fits.Column(name='flux', format='E', array=flux)]))
hdu.writeto('find_star_cands.fits')
'''

star_cand = fits.open('../find_star_cands.fits')[1].data
sgl = star_cand.gl.tolist()
sgb = star_cand.gb.tolist()
s_cgal = SkyCoord(sgl*u.degree, sgb*u.degree, frame='galactic')

###################################################################
# Load bstars
###################################################################
bstar = fits.open('../tychobstarmatch.fits')[1].data
bcut1 = np.where((bstar.nuv_cps > 10.) & (bstar.nuv_cps < 100.))
bstar1 = bstar[bcut1]
bcut2 = np.where((bstar1.gl > 0) & (bstar1.gl < 2) & (bstar1.gb > -2.) & (bstar1.gb < 0.))
bstar2 = bstar1[bcut2]
bgal2 = SkyCoord(bstar2.gl*u.degree, bstar2.gb*u.degree, frame='galactic')

###################################################################
# Match stars with detections
###################################################################
bstarind, star_candind, starangsep, stardist3d = search_around_sky(bgal2, s_cgal, 1.*u.arcmin)

uniqueind = np.unique(bstarind, return_index=True)[1]
anglim = np.where(starangsep.degree[uniqueind] < 0.002)

bstar3 = bstar2[bstarind[uniqueind]]
s_match = star_cand[star_candind[uniqueind]]

bstar3 = bstar3[anglim]
bstarname = bstar3.galex_id
bgl = bstar3.gl
bgb = bstar3.gb
bflux = bstar3.nuv_cps

s_match = s_match[anglim]
sgl = s_match.gl
sgb = s_match.gb
sflux = s_match.flux


# Make a table
table = 0
if table == 1:
	a1 = bstar3.galex_id
	a2 = bstar3.cat_id
	a3 = bstar3.ra
	a4 = bstar3.dec
	a5 = bgl
	a6 = bgb
	a7 = bstar3.fuv_cps
	a8 = bstar3.fuv_cat
	a9 = bstar3.nuv_cps
	a10 = bstar3.nuv_cat
	a11 = bstar3.BTmag
	a12 = bstar3.VTmag
	a13 = bstar3.VTmag - 0.090 * (bstar3.BTmag-bstar3.VTmag)
	a14 = a13 + 0.850 * (bstar3.BTmag-bstar3.VTmag)
	a15 = s_match.gl
	a16 = s_match.gb
	a17 = sflux
	col1 = fits.Column(name='galex_id', format='E', array=a1)
	col2 = fits.Column(name='cat_id', format='10A', array=a2)
	col3 = fits.Column(name='ra', format='E', array=a3)
	col4 = fits.Column(name='dec', format='E', array=a4)
	col5 = fits.Column(name='gl', format='E', array=a5)
	col6 = fits.Column(name='gb', format='E', array=a6)
	col7 = fits.Column(name='fuv_cps', format='E', array=a7)
	col8 = fits.Column(name='fuv_cat', format='10A', array=a8)
	col9 = fits.Column(name='nuv_cps', format='E', array=a9)
	col10 = fits.Column(name='nuv_cat', format='10A', array=a10)
	col11 = fits.Column(name='BTmag', format='E', array=a11)
	col12 = fits.Column(name='VTmag', format='E', array=a12)
	col14 = fits.Column(name='Bmag', format='E', array=a14)
	col13 = fits.Column(name='Vmag', format='E', array=a13)
	col15 = fits.Column(name='find_gl', format='E', array=a15)
	col16 = fits.Column(name='find_gb', format='E', array=a16)
	col17 = fits.Column(name='find_flux', format='E', array=a17)
	new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17])
	hdu = fits.BinTableHDU.from_columns(new_cols)
	hdu.writeto('../find_bstar_t2_matches.fits')
	

