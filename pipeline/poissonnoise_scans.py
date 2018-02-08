from astropy.io import fits, ascii
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
import photutils

scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '72-82', '81-91', '90-100', '9-19', '99-109']

imeds = []
cmeds = []
smoothed = 0

for scan in scans:
    print('Running '+scan)
    # Get cmed and imed
    hdu = fits.open('galexscans/count_map_'+scan+'_in.fits')
    img = hdu[0].data
    im1 = img[np.where(img > 0)]

    ct = fits.open('countscans/count_map_'+scan+'_count.fits')[0].data
    ct1 = ct[np.where(ct > 0)]
    
    cmed = np.median(ct1)
    imed = np.median(im1)
    del ct
    del ct1
    del im1
    
    imeds.append(imed)
    cmeds.append(cmed)
    
    '''
    medvals = Table.read('median_scanvals.txt', format='ascii')
    medvals = medvals[np.where(medvals['scans'] == scan)]
    imed = medvals['imed'][0]
    cmed = medvals['cmed'][0]

    # Find zero regions in img
    cutoffx, cutoffy = np.where(img < 10**-3)

    print('Now doing aperture photometry')
    apertures = photutils.CircularAperture((cutoffx, cutoffy), r=5)
    table = photutils.aperture_photometry(img, apertures)

    # Remove regions within data from list leaving only blank areas
    tcut = table[np.where(table['aperture_sum'] < 10**-3)]
    positions = (np.array(tcut['ycenter']).astype(int), np.array(tcut['xcenter']).astype(int))
    '''

    xmed = np.median(img, axis=0)

    if smoothed == 1:
        img = fits.open('smoothscans/im1_'+scan+'.fits')[0].data

    for line in range(len(img[0,:])):
        if xmed[line] < 10**-3:
            img[:,line] = imed*(1+(np.random.poisson(cmed, len(img[:,0])) - cmed)/cmed) 

    wcsmap = WCS(hdu[0].header)
    header = wcsmap.to_header()
    
    if smoothed == 0:
        fits.writeto('countmap_poisson/count_map_'+scan+'_in_poisson.fits', img, header)

    if smoothed == 1:
        fits.writeto('countmap_smoothed_poisson/im1_'+scan+'_poisson.fits', img, header)

maskvals = Table([imed, cmed], names=['imed', 'cmed'])
ascii.write(maskedvals, 'masked_values.txt', format='basic')
