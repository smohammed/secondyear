from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16,8
matplotlib.rcParams['font.size'] = 20

def plotmap(a,b,cutflux):
	a=a		# Start gl
	b=b		# End gl
	
	pixscale = 1 * (b - 20.)/20.
	
	#phomap = np.sqrt(fits.open('../intensitymap1200'+'_'+str(a)+'_'+str(b)+'.fits')[0].data).T
	phomap = np.sqrt(fits.open('../intmapcsvcorr1200'+'_'+str(a)+'_'+str(b)+'.fits')[0].data).T	
	#expmap = fits.open('../8-28-aspcorr_new_scst.fits')[0].data.T
	expmap = fits.open('../aspcorr_new_scst_1200.fits')[0].data.T

	bstars = fits.open('../bstar.fits')[1].data

	gl = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.l.degree
	gb = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.b.degree
	blim = np.where((bstars.fuv_cps > cutflux) & (gl > a) & (gl < b) & (gb > -10) & (gb < 10))

	gl = gl[blim]
	gb = gb[blim]

	offset = 0
	#expmap1 = expmap[39:1238,39+1199*pixscale:1238+1199*pixscale]
	expmap1 = expmap[39:1238,39+offset+1199*pixscale:1238+offset+1199*pixscale]
	
	intmap = phomap/expmap1	
	intmap = np.nan_to_num(intmap)

	plt.scatter(gl,gb,facecolor='none',edgecolor='red')
	plt.xlabel('gl')
	plt.ylabel('gb')
	return plt.imshow(intmap,vmin=0,vmax=0.4,origin='lower',extent=[a,b,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray), plt.show()	

