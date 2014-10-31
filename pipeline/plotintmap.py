from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np

def plotmap(a,b):
	a=a		# Start gl
	b=b		# End gl
	
	pixscale = 1 * (b - 20.)/20.
	
	phomap = np.sqrt(fits.open('../intensitymap1200'+'_'+str(a)+'_'+str(b)+'.fits')[0].data).T
	expmap = fits.open('../8-28-aspcorr_new_scst.fits')[0].data.T
	
	expmap1 = expmap[39:1238,39+1199*pixscale:1238+1199*pixscale]
	
	intmap = phomap/expmap1
	
	return plt.imshow(intmap,vmin=0,vmax=0.4,origin='lower',extent=[a,b,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray), plt.show()	

