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
#scans = ['0-10', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109']

scans = ['144-154']

skyrange = scans

weight = 0

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
    if weight == 0:
        os.system('sextractor ../../galexscans/im1_'+region+'_in_03_20_18.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im1_'+region+'.fits')
    
    if weight == 1:
        os.system('sextractor ../../galexscans/im1_'+region+'_in_03_20_18.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im1_'+region+'_weighted.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ../../galexscans/count_map_'+region+'_exp_noise_04-18-18.fits')

    print 'sextractor run 1 finished'

    if weight == 0:
        bkgd = fits.open('../../galexscans/background_im1_'+region+'.fits')[0].data

    if weight == 1:
        bkgd = fits.open('../../galexscans/background_im1_'+region+'_weighted.fits')[0].data

    im1 = img - bkgd
    header = wcsmap.to_header()

    try:
        fits.writeto('../../galexscans/im1_'+region+'_bksub.fits', im1, header, overwrite=True)

    except IOError:
        os.remove('../../galexscans/im1_'+region+'_bksub.fits')
        fits.writeto('../../galexscans/im1_'+region+'_bksub.fits', im1, header, overwrite=True)

    # With no background step, subtract background prior to this
    if weight == 0:
        os.system('sextractor ../../galexscans/im1_'+region+'_bksub.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0')

    if weight == 1:
        os.system('sextractor ../../galexscans/im1_'+region+'_bksub.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ../../galexscans/count_map_'+region+'_exp_noise_04-18-18.fits')

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

    scatter_contour(alldata['nuv'], alldata['FWHM_IMAGE'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10))
    plt.xlim(12, 21)
    plt.xlabel('NUV')
    plt.ylabel('FWHM [pixels]')
    plt.title('')
    plt.savefig('fwhmvsnuv_tests_detect_thresh_6.png')

    if weight == 0:
        ascii.write(alldata, '../../galexscans/tests_starcat_'+currregion+'_11_08_2018.txt', format='ipac')

    if weight == 1:
        ascii.write(alldata, '../../galexscans/starcat_'+currregion+'_weighted_04_18_2018.txt', format='ipac')


    print 'Converted to gl, gb, finished'
