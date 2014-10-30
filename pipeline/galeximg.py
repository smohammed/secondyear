from astropy.io import fits
import numpy as np
from numpy import loadtxt
#from matplotlib import pyplot as plt

########################################################################################
# This code generates an intensity map. Combine with the exposure map to get a final image.
########################################################################################

############################################
# Input galactic longitude range
############################################
a=120 		# Start gl
b=140		# End gl

############################################
# Set the file pathways
############################################

#path = "http://uqbar.astro.columbia.edu/~bsktam/galactic_plane/fits2/"
path = "../intmapfiles/"
lines = loadtxt("../filelists/fits_photon_list.txt", comments = "#", delimiter = ",", unpack = False, dtype = str)

############################################
# Grab appropriate files using gl range
############################################

startgl = a
endgl = a

while (float(lines[startgl].split('_')[0])/10.) < a:
	startgl+=1
while (float(lines[endgl].split('_')[0])/10.) < b:
	endgl+=1
endgl = endgl - 1


photons = fits.open(path+lines[startgl])[1].data
photons2 = fits.open(path+lines[startgl+1])[1].data

gl = np.concatenate((photons.gl_cor,photons2.gl_cor))
gb = np.concatenate((photons.gb_cor,photons2.gb_cor))

############################################
# Combine photons from up and down runs, stack data horizontally
############################################

for i in range(startgl+2,endgl,2):
	print i
	try:
		f1 = fits.open(path+lines[i])[1].data
		f1gl_cor = f1.gl_cor
		f1gb_cor = f1.gb_cor

	except IOError:								# if the file is corrupted, skip it 
		print 'Index '+str(i)+' is corrupted'

	try:
		f2 = fits.open(path+lines[i+1])[1].data
		f2gl_cor = f2.gl_cor
		f2gb_cor = f2.gb_cor

	except IOError:								# if the file is corrupted, skip it 
		print 'Index '+str(i+1)+' is corrupted'
	
	glinit = np.concatenate((f1gl_cor,f2gl_cor))
	gbinit = np.concatenate((f1gb_cor,f2gb_cor))
	gl = np.hstack([gl, glinit])
	gb = np.hstack([gb, gbinit])

############################################
# Bin data and plot stuff
############################################

binnum = 1200

H, xbins, ybins = np.histogram2d(gl, gb, bins = (np.linspace(a, b, binnum), np.linspace(-10, 10, binnum)))
#fig = plt.figure(figsize = (10,10))
#ax = plt.axes()

fits.PrimaryHDU(H).writeto('../intensitymap'+str(binnum)+'_'+str(a)+'_'+str(b)+'.fits')

'''
ax.imshow(np.sqrt(H).T, vmin = 0, vmax = 5, origin = 'lower', extent = [a, b, -10, 10], interpolation = 'nearest', aspect = 'auto', cmap = 'gray')
ax.set_xlabel(r'${\rm gl}$')
ax.set_ylabel(r'${\rm gb}$')
ax.set_title('Galactic Plane in UV, %d<gl<%d'%(a,b))

fig.set_dpi(300)
fig.savefig('GALEX_galplane_%d-%d_vmax05.png'%(a,b), bbox_inches = 'tight', pad_inches = 0)
'''