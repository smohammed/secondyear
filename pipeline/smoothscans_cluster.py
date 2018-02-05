from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import os
from astropy.convolution import convolve, Gaussian2DKernel

# Files written:
# 1. Cutouts of main image to remove edges
# 2. SExtractor output, background used
# 3. starcatalog*.fits which adds WCS data to stars

#########################################################################
# Select desired field from list
#########################################################################
# Already run
#scans =['0-10', '18-28', '270-280', '351-1', '108-118', '189-199', '27-37', '36-46', '117-127', '198-208', '279-289', '45-55', '126-136', '207-217', '288-298', '63-73', '135-145', '297-307', '72-82', '144-154', '225-235', '306-316', '81-91', '153-163', '234-244', '315-325', '90-100', '162-172', '243-253', '324-334', '9-19', '171-181', '252-262', '333-343', '99-109', '180-190', '261-271', '342-352', '216-226']


scans = ['225-235']

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
    print(currregion)
    region = currregion.replace('.', '')

    hdu = fits.open('countscans/count_map_'+region+'_count.fits')[0]
    img = hdu.data
    wcsmap = WCS(hdu.header)
    '''
    #find boundry
    exp_hdu_list = fits.open('galexscans/count_map_'+region+'_exp.fits')
    exp = exp_hdu_list[0].data
    exp[np.isnan(exp)] = 0.
    exp[exp>1000] = 0.
    exp[exp<0] = 0.
    ver = np.sum(exp, axis=0)
    hor = np.sum(exp, axis=1)
    half = np.argmax(ver)
    sort =  np.argsort(np.absolute(ver-np.max(ver)*0.2))
    im1xmin = sort[sort<half][0]
    im1xmax = sort[sort>half][0]
    print im1xmin, im1xmax, np.absolute(im1xmin-im1xmax)
    half = np.argmax(hor)
    sort =  np.argsort(np.absolute(hor-np.max(hor)*0.2))
    im1ymin = sort[sort<half][0]
    im1ymax = sort[sort>half][0]
    print im1ymin, im1ymax, np.absolute(im1ymin-im1ymax)
    exp_hdu_list.close()
    '''

    #########################################################################
    # Make cutouts of initial image to help with background correction
    #########################################################################
    #img = img[im1ymin:im1ymax, im1xmin:im1xmax]
    #wcsmap = wcsmap[im1ymin:im1ymax, im1xmin:im1xmax]
    header = wcsmap.to_header()

    gauss = Gaussian2DKernel(stddev=3)
    im1 = convolve(img, gauss)
    print('Smoothing finished')

    try:
        fits.writeto('smoothcounts/im1_'+region+'_counts.fits', im1, header, clobber=True)

    except IOError:
        os.remove('smoothcounts/im1_'+region+'_counts.fits')
        fits.writeto('smoothcounts/im1_'+region+'_counts.fits', im1, header, clobber=True)

