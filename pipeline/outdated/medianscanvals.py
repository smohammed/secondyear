from astropy.io import fits, ascii
import numpy as np

scans = ['63-73', '180-190', '225-235']
imed = []
cmed = []

for scan in scans:
	ct = fits.open('im1_'+scan+'_count.fits')[0].data
	hdu = fits.open('im1_'+scan+'.fits')
	img = hdu[0].data

	cmed1 = np.median(ct[np.where(ct > 0)])
	imed1 = np.median(img[np.where(img > 10**-5)])

	imed.append(imed1)
	cmed.append(cmed1)


table = Table([scans, imed, cmed], names=['scans', 'imed', 'cmed'])	
ascii.write(table, 'median_scanvals.txt', format='basic')

