from astropy.io import fits, ascii
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
import scipy.ndimage as ndimage

# Do these later
scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '72-82', '81-91', '90-100', '9-19', '99-109']

smoothed = 1

for scan in scans:
    print('Running '+scan)
    # Gather cmean and imean
    hdu = fits.open('galexscans/count_map_'+scan+'_in.fits')
    img = hdu[0].data

    if smoothed == 0:
        meanvals = Table.read('mean_scanvals.txt', format='ascii')

    if smoothed == 1:
        meanvals = Table.read('mean_scanvals_smoothed.txt', format='ascii')

    meanvals = meanvals[np.where(meanvals['scan'] == scan)]
    imean = meanvals['imean'][0]
    cmean = meanvals['cmean'][0]

    if smoothed == 1:
        del img
        img = fits.open('smoothscans/im1_'+scan+'.fits')[0].data

    nval = imean * (np.random.poisson(cmean, 100000)/cmean)
    struc3 = ndimage.generate_binary_structure(2,2).astype(img.dtype)
    mask = ~ndimage.binary_dilation(img, structure=struc3, iterations=1)
    np.place(img, mask, nval.astype(np.float32))

    wcsmap = WCS(hdu[0].header)
    header = wcsmap.to_header()

    if smoothed == 0:
        fits.writeto('countmap_poisson/count_map_'+scan+'_in_poisson_03_13_18.fits', img, header)

    if smoothed == 1:
        fits.writeto('countmap_smoothed_poisson/im1_'+scan+'_poisson_03_13_18.fits', img, header, ove
rwrite=True)
