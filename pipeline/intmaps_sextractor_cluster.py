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
#scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '189-199', '198-208', '207-217', '216-226', '270-280', '279-289', '288-298', '297-307', '351-1']

# Run sextractor on these
scans = ['144-154', '153-163', '162-172', '171-181', '180-190', '225-235', '234-244']

# Scans that still need smoothing
#scans = ['243-253', '252-262', '261-271', '306-316', '315-325', '324-334', '333-343', '342-352', '72-82', '81-91', '90-100', '9-19', '99-109']

skyrange = scans

# Custom parameters that I use
#filt = 'mexhat_4.0_9x9'
filt = 'gauss_3.0_7x7.conv'
code = 'det_thresh4_phot_autopar2.5_3.5'
# Det Thresh = 4, Analysis Thresh = 3.5

'''
def createCircularMask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask
'''

#########################################################################
# Decide which run to do
#########################################################################
run1 = 1
run2 = 0

#########################################################################
# Choose a field as defined above
#########################################################################
for currregion in skyrange:
    #########################################################################
    # Figure out what the cutout range should be
    #########################################################################
    print 'current region = ' + currregion
    region = currregion.replace('.', '')

    hdu = fits.open('../../galexscans/count_map_'+region+'_in.fits')[0]
    img = hdu.data
    wcsmap = WCS(hdu.header)
    '''
    #find boundry
    exp_hdu_list = fits.open('../../galexscans/count_map_'+region+'_exp.fits')
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
    # Run sextractor, subtract background from original and run again
    #########################################################################
    os.system('sextractor ../../galexscans/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im1_'+region+'.fits')

    bkgd = fits.open('../../galexscans/background_im1_'+region+'.fits')[0].data
    im1 = img - bkgd
    header = wcsmap.to_header()     

    try:
        fits.writeto('../../galexscans/im1_'+region+'.fits', im1, header, clobber=True)

    except IOError:
        os.remove('../../galexscans/im1_'+region+'.fits')
        fits.writeto('../../galexscans/im1_'+region+'.fits', im1, header, clobber=True)

    # With no background step, subtract background prior to this
    os.system('sextractor ../../galexscans/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0')

    print 'SExtractor finished'

    #########################################################################
    # Get output from sextractor, convert NUV
    #########################################################################
    im1sex = Table.read('../../galexscans/sex_im1_'+region+'.fits', format='fits')

    data = im1sex
    nuv = -2.5*np.log10(data['FLUX_AUTO']) + 20.08
    data['nuv'] = nuv
    data = data[~np.isnan(data['nuv'])]

    #########################################################################
    # Fix if first scan, save table
    #########################################################################
    skygal = SkyCoord(data['ALPHA_J2000'], data['DELTA_J2000'], frame='icrs').galactic
    coord = Table([skygal.l.degree, skygal.b.degree], names=('gl', 'gb'))
    if region == '5':
        coord['gl'][np.where(coord['gl'] > 350)] = coord['gl'][np.where(coord['gl'] > 350)] - 360

    alldata = hstack([data, coord])

    ascii.write(alldata, '../../galexscans/starcat_'+currregion+'_01_10_2018.txt', format='ipac')

    print 'Converted to gl, gb, finished'
