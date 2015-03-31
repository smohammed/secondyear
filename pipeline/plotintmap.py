from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16,8
matplotlib.rcParams['font.size'] = 20

def plotmap(startgl,endgl,startgb,endgb,cutfluxlow, cutfluxhigh):
	a = startgl		# Start gl
	b = endgl		# End gl
	
	c = startgb
	d = endgb

	###################################################################
	# Set a scale so we can cut the correct portion of the sky from the big expmap 
	###################################################################
	pixscalex = (b - 20.)/20.
	pixscaley = (d)

	offset = 0

	###################################################################
	# Load the photon map, which is 20x20 deg
	###################################################################
	#phomap = np.sqrt(fits.open('../intensitymap1200'+'_'+str(a)+'_'+str(b)+'.fits')[0].data).T
	phomap = np.sqrt(fits.open('../intmapcsvcorr1200'+'_'+str(a)+'_'+str(b)+'.fits')[0].data).T	

	###################################################################
	# Load the expmap which is 360x20 deg
	###################################################################
	expmap = fits.open('../12-18-aspcorr_new_scst_1200.fits')[0].data.T

	###################################################################
	# Set pixel limits on the image to remove the edges of the expmap, then apply to the map
	###################################################################
	ymin = 39 # is gb = -10 
	ymax = 1238 # is gb = 10 
	xmin = 39 + 1199 * pixscalex + offset 
	xmax = 1238 + 1199 * pixscalex + offset

	#expmap1 = expmap[39:1238,39+1199*pixscalex:1238+1199*pixscalex]
	#expmap1 = expmap[39:1238,39+offset+1199*pixscalex:1238+offset+1199*pixscalex]
	expmap1 = expmap[ymin:ymax,xmin:xmax]	

	###################################################################
	# Find gl,gb of bstars in that region, restrict by flux
	###################################################################
	bstars = fits.open('../bstar.fits')[1].data
	gl = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.l.degree
	gb = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.b.degree
	blim = np.where((bstars.nuv_cps > cutfluxlow) & (bstars.nuv_cps < cutfluxhigh) & (gl > a) & (gl < b) & (gb > -10) & (gb < 10))
	gl = gl[blim]
	gb = gb[blim]

	###################################################################
	# Combine the exposure map with the intensity map
	###################################################################
	intmap = phomap/expmap1	
	intmap = np.nan_to_num(intmap)

	###################################################################
	# Now plot!
	###################################################################
	plt.scatter(gl,gb,facecolor='none',edgecolor='red')
	plt.xlabel('gl')
	plt.ylabel('gb')
	return plt.imshow(intmap,vmin=0,vmax=0.7,origin='lower',extent=[a,b,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray), plt.show()	

