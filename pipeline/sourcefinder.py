# Matches b stars with GALEX photons

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky

##############################################################
##############################################################
# TO DO
# 1. Fix stupid table error with times and nphoton
# 2. Automate for all slices. Remove duplicate lines and create lists for things with similar times
# 3. Find the average shift of photons from a star, that will indicate the needed offset.
##############################################################
##############################################################

##############################################################
# Load B stars and make NUV cut for 10-100 cts/s (for 20-15 mag)
##############################################################
bstar = fits.open('../bstar.fits')[1].data
bcut1 = np.where((bstar.nuv_cps > 10.) & (bstar.nuv_cps < 100.))
bstar = bstar[bcut1]

##############################################################
# Load CSV files to match up times (cut to one file for now), fix gl values
##############################################################
path = "../corrcsv/"
lines = np.loadtxt("../filelists/csv_galcorr.txt", dtype = 'string')

field = fits.open('../corrcsv/00005_0001galcorr.csv.fits')[1].data

#for skyfield in range(len(lines)):
#	field = fits.open(path+lines[skyfield])

##############################################################
# Cut into one second slices & limit bstar data range
##############################################################
timestep = np.arange(np.round(field.time[-1])+1)

names,times,nphotons,stargl,stargb,phogl,phogb = [],[],[],[],[],[],[]

# Limit to the first three timesteps to make it run faster
for tstep in timestep[:3]:
	fieldlim = np.where((field.time >= tstep) & (field.time < tstep+1.))
	fgl = field.gl[fieldlim]
	fgb = field.gb[fieldlim]

	# Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
	for i in range(len(fgl)):
		if min(fgl) < 10. and fgl[i] > 310.:
			fgl[i] = fgl[i] - 360.
	
	# Convert to galactic coordinates and limit the bstar range 
	fgal = SkyCoord(fgl.tolist()*u.degree, fgb.tolist()*u.degree, frame='galactic')
	bgal = SkyCoord(bstar.ra*u.degree, bstar.dec*u.degree, frame='icrs').galactic
	bcut2 = np.where((bgal.l.degree > np.min(fgl)) & (bgal.l.degree < np.max(fgl)) & (bgal.b.degree > -10.) & (bgal.b.degree < 10.))
	bstar2 = bstar[bcut2]
	bgal = bgal[bcut2]

	##############################################################
	# Find photons within 2' of each star 
	##############################################################
	starind, photind, starsep, stardist = search_around_sky(bgal,fgal,2.*u.arcmin)
	starname = bstar2.galex_id[starind]
	bgl = bgal.l.degree[starind]
	bgb = bgal.b.degree[starind]
	
	# Do the same gl fix as for the fgl values
	for i in range(len(bgl)):
		if min(bgl) < 10. and bgl[i] > 310.:
			bgl[i] = bgl[i] - 360.
	
	# Check to see that you're not totally screwing up by plotting
	plot = 0
	if plot == 1:
		plt.scatter(fgl[photind],fgb[photind],c='orange',s=10,edgecolor='none')
		plt.scatter(bgl,bgb,c='green',s=20)
		#plt.scatter(fgl,fgb,c='red',s=3,edgecolor='none')
		#plt.show()
		#break

	# Add photon times for each index
	photon = []
	for i in starind:
		photon.append(sum(starind == i))


	




	#Add values to arrays so you can then make a big table
	cutnames = np.unique(starname)
	bgl = np.unique(bgl)
	bgb = np.unique(bgb)
	fieldgl = fgl[photind].tolist()
	fieldgb = fgb[photind].tolist()


	print cutnames
	print fieldgl
	print len(starname)
	print len(fieldgl)
	print np.average(fieldgl)
	print np.average(fieldgb)

	names = np.concatenate((names,starname),axis=0)
	stargl = np.concatenate((stargl,bgl),axis=0)
	stargb = np.concatenate((stargb,bgb),axis=0)
	phogl.append(fieldgl)
	phogb.append(fieldgb)
	times = np.concatenate((times,np.ones(len(names))*tstep),axis=0)
	nphotons = np.concatenate((nphotons,photon),axis=0)


##############################################################
# Now find average photon displacement from each star
##############################################################





##############################################################
# Make the table
##############################################################
make_table = 0
if make_table == 1:
	a1 = np.array(names)
	a2 = np.array(stargl)
	a3 = np.array(stargb)
	a4 = np.array(phogl)
	a5 = np.array(phogb)
	a6 = np.array(times)
	a7 = np.array(nphotons)
	cols = []
	col1 = fits.Column(name='name', format='E', array=a1)
	col2 = fits.Column(name='star_gl', format='E', array=a2)
	col3 = fits.Column(name='star_gb', format='E', array=a3)
	col4 = fits.Column(name='photon_gl', format='20E', array=a4)
	col5 = fits.Column(name='photon_gb', format='E', array=a5)
	col6 = fits.Column(name='times', format='E', array=a6)
	col7 = fits.Column(name='nphotons', format='E', array=a7)
	cols.append(col1)
	cols.append(col2)
	cols.append(col3)
	cols.append(col4)
	cols.append(col5)
	cols.append(col6)
	cols.append(col7)
	new_cols = fits.ColDefs(cols)
	hdu = fits.BinTableHDU.from_columns(new_cols)
	hdu.writeto('00005_0001galcorr_bstarphotons.fits')


