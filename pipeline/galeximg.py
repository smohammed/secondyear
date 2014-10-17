from astropy.io import fits
import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

########################################################################################
# This code generates an intensity map. Combine with the exposure map to get a final image.
########################################################################################

############################################
# Input parameters
############################################
a=0 		# Start gl
b=20		# End gl
x=1 		# 
y=612 		# 

############################################
# Set the file pathways
############################################

#path = "http://uqbar.astro.columbia.edu/~bsktam/galactic_plane/fits2/"
path = "../intmapfiles/"
lines = loadtxt("../filelists/fits_photon_list.txt", comments = "#", delimiter = ",", unpack = False, dtype = str)

############################################
# Probably not important
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
# 
############################################

for i in range(startgl+2,endgl,2):
	print i
	try:
		f1 = fits.open(path+lines[i])[1].data
		f2 = fits.open(path+lines[i+1])[1].data
		f1gl_cor = f1.gl_cor
		f1gb_cor = f1.gb_cor
		f2gl_cor = f2.gl_cor
		f2gb_cor = f2.gb_cor[::-1]

	except IOError:
		print 'Index '+str(i)+' is corrupted'
		#f1gl_cor = np.zeros(np.size())

	#NOPE figure out how to add the arrays first correctly!	

	glinit = np.concatenate((f1gl_cor,f2gl_cor))
	gbinit = np.concatenate((f1gb_cor,f2gb_cor))
	gl = np.hstack([gl, glinit])
	gb = np.hstack([gb, gbinit])


############################################
# Plot stuff
############################################

H, xbins, ybins = np.histogram2d(gl, gb, bins = (np.linspace(a, b, 1200), np.linspace(-10, 10, 1200)))
fig = plt.figure(figsize = (10,10))
ax = plt.axes()

#fits.PrimaryHDU(H).writeto('intensitymap.fits')

ax.imshow(np.sqrt(H).T, vmin = 0, vmax = 5, origin = 'lower', extent = [a, b, -10, 10], interpolation = 'nearest', aspect = 'auto', cmap = 'gray')
ax.set_xlabel(r'${\rm gl}$')
ax.set_ylabel(r'${\rm gb}$')
ax.set_title('Galactic Plane in UV, %d<gl<%d'%(a,b))

fig.set_dpi(300)
fig.savefig('GALEX_galplane_%d-%d_vmax05.png'%(a,b), bbox_inches = 'tight', pad_inches = 0)
