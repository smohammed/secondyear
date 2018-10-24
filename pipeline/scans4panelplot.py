from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import matplotlib, matplotlib.cm as cm
from matplotlib import pyplot as plt
matplotlib.rcParams['figure.figsize'] = 12, 12

scans = ['0-10', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109']

for scan in scans:
	print 'scan = '+ scan
	hdu = fits.open('../../galexscans/count_map_'+scan+'_in_poisson_03_20_18.fits')[0]
	wcs = WCS(hdu.header)
	bkgd = fits.open('../../galexscans/background_im1_'+scan+'.fits')[0].data
	t = Table.read('../../galexscans/starcat_'+scan+'_03_26_2018.txt', format='ascii')

	fig = plt.figure()
	fig.add_subplot(221, projection = wcs)
	plt.imshow(hdu.data, origin='lower', cmap=cm.gray, vmin=0, vmax=0.1, aspect='auto')

	fig.add_subplot(222, projection=wcs)
	plt.imshow(bkgd, origin='lower', cmap=cm.gray, vmin=0, vmax=0.1, aspect='auto')

	fig.add_subplot(223)
	plt.scatter(t['gl'], t['gb'], c=t['FWHM_IMAGE'], vmin=0, vmax=10, s=3)
	plt.colorbar().set_label('FWHM')
	plt.gca().invert_xaxis()

	fig.add_subplot(224)
	plt.scatter(t['nuv'], t['FWHM_IMAGE'], s=1, c='black', marker='.')
	plt.ylim(0, 10)
	plt.xlabel('NUV')
	plt.ylabel('FWHM')
	plt.savefig('03-26-18_scan4plots_'+scan+'.png')
	plt.clf()
	plt.close()
