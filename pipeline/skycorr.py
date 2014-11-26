from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

csvpath = "../csv/"
csvfiles = np.loadtxt("../filelists/AIS_GAL_SCAN_csv.txt", comments = "#", delimiter = ",", unpack = False, dtype = str)

scstpath = "../scst/"
scstfiles = np.loadtxt('../filelists/aspcorr_new_list.txt',dtype='string')


for csvline in range(len(csvfiles)):
	# Load photon file
	pho = fits.open(csvpath+csvfiles[csvline])[1].data

	# Find correct scst file for the photon file
	for scstline in range(len(scstfiles)):						# Find correct photon times file
		if csvfiles[csvline][13:23] == scstfiles[scstline][13:23]:
			scstnum = scstline

	# Open the scst file and define needed variables
	exp = fits.open(scstpath+scstfiles[scstnum])[1].data
	scsttime = exp['T']
	dglcrefine = exp.DGL_CREFINE
	dgbcrefine = exp.DGB_CREFINE

	# Set time limit for data to remove useless junk
	timelim = np.where(pho['T'] > 0)
	pho = pho[timelim]

	# normalize times so they start at 0 and not at weird numbers
	photime = (pho['T'] - pho['T'][0])/1000.
	exptime = scsttime - scsttime[0]

	# convert to galactic coordinates
	gl = SkyCoord(pho.RA*u.degree, pho.DEC*u.degree, frame='icrs').galactic.l.degree
	gb = SkyCoord(pho.RA*u.degree, pho.DEC*u.degree, frame='icrs').galactic.b.degree


	# Match photon times with scst times, apply gl/gb correction
	for i in range(len(exptime)-1):
		cuts = (photime >= exptime[i]) & (photime < exptime[i+1])
		gl[cuts] = gl[cuts] - exp.DGL_CREFINE[i]
		gb[cuts] = gb[cuts] - exp.DGB_CREFINE[i]

	# Because I can't write code well, fix last value
	cuts = (photime == exptime[-1])
	gl[cuts] = gl[cuts] - exp.DGL_CREFINE[-1]
	gb[cuts] = gb[cuts] - exp.DGB_CREFINE[-1]

	# MAKE THE ARRAY
	a1 = np.array(photime)
	a2 = np.array(gl)
	a3 = np.array(gb)
	cols = []
	col1 = fits.Column(name='time', format='E', array=a1)
	col2 = fits.Column(name='gl', format='E', array=a2)
	col3 = fits.Column(name='gb', format='E', array=a3)
	cols.append(col1)
	cols.append(col2)
	cols.append(col3)
	
	#orig_cols = pho.columns
	new_cols = fits.ColDefs(cols)
	#hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
	hdu = fits.BinTableHDU.from_columns(new_cols)
	hdu.writeto('../corrcsv/'+csvfiles[csvline][13:23]+'galcorr.csv.fits')	