from astropy.io import fits
import numpy as np

# These can be arbitrary 

nx=81.
ny=81.

#################################################
# Make array exposure map
#################################################
# add a buffer around edge just for good measure
#expim=np.zeros((360.0*60+nx,20.0*60+ny))

expim=np.zeros((20.0*60+nx,20.0*60+ny)) # For 20x20 deg field

#expim=np.zeros((200.0*60+nx,200.0*60+ny))

#################################################
# Creates a circular mask array of 1s bordered by 0s
#################################################
mask = np.zeros((nx+1.,ny+1.))
radius = (1.24/2)*60.

for i in np.arange(nx):
     for j in np.arange(ny):
          a = i-40.
          b = j-40.
          if np.sqrt(a**2+b**2) < radius:
               mask[i][j] = 1.

#################################################
# Open files to mask
#################################################

#files1 = np.loadtxt('../filelists/aspcorr_new_list.txt',dtype='string')
files1 = np.loadtxt('../filelists/truncaspcorr_new_list.txt',dtype='string') # For 0-20deg
path = "../scst/"
photonlist = np.loadtxt('../filelists/photontimes.txt',dtype='string')


for line in range(len(files1)):									# Iterate through scst files
	print line
	d = fits.open(path+files1[line])[1].data					# Open exp file
	for pholine in range(len(photonlist)):						# Find correct photon times file
		if files1[line].split('_')[3] == photonlist[pholine][:5]:
			pfile = pholine

	photimes = np.loadtxt('../times/'+photonlist[pfile],comments='#',delimiter=',',unpack=False,dtype=str)[:-1]													# Load correct photon times file

	delpho=[]
	for i in range(len(photimes)-1):
 	   delpho.append(np.where((d['T']-d['T'][0])==float(photimes[i])))

 	q=[]
 	for i in range(len(delpho)):
 		if len(delpho[i][0]) > 0:
 			q.append(delpho[i][0][0])


	#################################################
	# Convert RA, DEC to gl, gb then to pixels 
	#################################################
	
	##### CONVERT RA, DEC to GL, GB, and then to PIXELS, assuming that each pixel is 1 arcmin (thus the x60 below)	
	##### THESE ARE ARRAYS of SCAN position at a given time (in 1 second intervals)
	
	#gl = SkyCoord(d.ra_acs*u.degree, d.dec_acs*u.degree, frame='icrs').galactic.l.degree
	#gb = SkyCoord(d.ra_acs*u.degree, d.dec_acs*u.degree, frame='icrs').galactic.b.degree
	
	gl = d.SCAN_GL[q]
	gb = d.SCAN_GB[q]

	gx = gl * 60.
	gy = (gb + 10.) * 60.

		
	for i in range(len(gx)):
		if gx[i] < 0:
			gx[i] = 0
		elif gx[i] > 360*60.:
			gx[i] = 360*60.

	for i in range(len(gy)):
		if gy[i] < 0:
			gy[i] = 0
		elif gy[i] > 20*60.:
			gy[i] = 20*60.

	##### ADD BORDER TO PIXEL VALUE
	# use 40 offset to keep from trailing off to account for errors at the border of the map
	gx = gx + 40.
	gy = gy + 40.
	
	#################################################
	# Apply mask and only take data with > 15000 photon counts
	#################################################
	
	##### ONLY USE SCAN POSITIONS FROM TIMES WHERE WE RECEIVED MORE THAN 15000 photons per second.  You can do something similar by only taking scan positions that have
	##### photons that you used to build the images.
	
	for j in range(len(d[q])-1):
		if d.NDCTEC[q][j] > 15000.:
			if np.shape(expim[gx[j]-41.:gx[j]+41., gy[j]-41.:gy[j]+41.]) == np.shape(mask):
				expim[gx[j]-41.:gx[j]+41., gy[j]-41.:gy[j]+41.] = expim[gx[j]-41.:gx[j]+41., gy[j]-41.:gy[j]+41.] + mask 
			else:
				print 'shape of expim != shape of mask' 

hdu = fits.PrimaryHDU(expim)
hdu.writeto('../12-18-aspcorr_new_scst_1200_0to20.fits')