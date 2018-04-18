from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import os

# Files written:
# 1. Background 
# 2. Background subtracted image to run sextractor in step 2
# 3. SExtractor output
# 4. starcatalog*.fits which adds WCS data to stars

#########################################################################
# Select desired field from list
#########################################################################
# Total scans
# To do
#scans = ['0-10', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', 

scans = ['234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109']

#scans = ['63-73']

skyrange = scans

# Custom parameters that I use
#filt = 'mexhat_4.0_9x9'
filt = 'gauss_3.0_7x7.conv'
code = 'det_thresh4_phot_autopar2.5_3.5'
# Det Thresh = 4, Analysis Thresh = 3.5

#########################################################################
# Choose a field as defined above
#########################################################################
for currregion in skyrange:
    #########################################################################
    # Figure out what the cutout range should be
    #########################################################################
    print 'current region = ' + currregion
    region = currregion.replace('.', '')

    hdu = fits.open('../../galexscans/count_map_'+region+'_in_poisson_03_20_18.fits')[0]
    img = hdu.data
    wcsmap = WCS(hdu.header)

    #########################################################################
    # Run sextractor, subtract background from original and run again
    #########################################################################
    os.system('sextractor ../../galexscans/im1_'+region+'_in_03_20_18.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im1_'+region+'.fits')

    print 'sextractor run 1 finished'

    bkgd = fits.open('../../galexscans/background_im1_'+region+'.fits')[0].data
    im1 = img - bkgd
    header = wcsmap.to_header()

    try:
        fits.writeto('../../galexscans/im1_'+region+'_bksub.fits', im1, header, overwrite=True)

    except IOError:
        os.remove('../../galexscans/im1_'+region+'_bksub.fits')
        fits.writeto('../../galexscans/im1_'+region+'_bksub.fits', im1, header, overwrite=True)

    # With no background step, subtract background prior to this
    os.system('sextractor ../../galexscans/im1_'+region+'_bksub.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0')

    print 'SExtractor run 2 finished'

    #########################################################################
    # Get output from sextractor, convert NUV
    #########################################################################
    im1sex = Table.read('../../galexscans/sex_im1_'+region+'.fits', format='fits')

    data = im1sex
    nuv = -2.5*np.log10(data['FLUX_AUTO']) + 20.08
    data['nuv'] = nuv
    data = data[~np.isnan(data['nuv'])]
    data = data[np.where(data['FWHM_IMAGE'] > 0)]

    #########################################################################
    # Fix if first scan, save table
    #########################################################################
    skygal = SkyCoord(data['ALPHA_J2000'], data['DELTA_J2000'], frame='icrs').galactic
    coord = Table([skygal.l.degree, skygal.b.degree], names=('gl', 'gb'))
    if region == '5':
        coord['gl'][np.where(coord['gl'] > 350)] = coord['gl'][np.where(coord['gl'] > 350)] - 360

    alldata = hstack([data, coord])

    ascii.write(alldata, '../../galexscans/starcat_'+currregion+'_03_26_2018.txt', format='ipac')

    print 'Converted to gl, gb, finished'
