from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
import photutils

scans = ['180-190', '225-235', '63-73']

for scan in scans:
    print('doing '+scan)
    # Get cmed and imed
    hdu = fits.open('smoothscans/im1_'+scan+'.fits')
    img = hdu[0].data

    medvals = Table.read('median_scanvals.txt', format='ascii')
    medvals = medvals[np.where(medvals['scans'] == scan)]
    imed = medvals['imed'][0]
    cmed = medvals['cmed'][0]

    n = imed * (1+(np.random.poisson(cmed*10, img.shape)/10. - cmed)/cmed)

    # Find zero regions in img
    cutoffx, cutoffy = np.where(img < 10**-3)

    apertures = photutils.CircularAperture((cutoffx, cutoffy), r=5)
    table = photutils.aperture_photometry(img, apertures)

    # Remove regions within data from list leaving only blank areas
    tcut = table[np.where(table['aperture_sum'] < 10**-3)]
    positions = (np.array(tcut['ycenter']).astype(int), np.array(tablecut['xcenter']).astype(int))

    # Multiply zeromask by n then finally add to img
    zeromask = np.zeros(img.shape)
    zeromask[positions] = True
    z = zeromask * n
    im2 = img + z

    wcsmap = WCS(hdu[0].header)
    header = wcsmap.to_header()
    fits.writeto('im1_'+scan+'_maskedpoisson.fits', im2, header)
