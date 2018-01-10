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
#scans = ['216-226']

scans =['0-10', '18-28', '270-280', '351-1', '108-118', '189-199', '27-37', '36-46', '117-127', '198-208', '279-289', '45-55', '126-136', '207-217', '288-298', '63-73', '135-145', '297-307']

# Run these later
#'72-82', '144-154', '225-235', '306-316', '81-91', '153-163', '234-244', '315-325', '90-100', '162-172', '243-253', '324-334', '9-19', '171-181', '252-262', '333-343', '99-109', '180-190', '261-271', '342-352']

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
# Decide which run
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

    hdu = fits.open('galexscans/count_map_'+region+'_in.fits')[0]
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
    if run1 == 1:
        img = img[im1ymin:im1ymax, im1xmin:im1xmax]
        wcsmap = wcsmap[im1ymin:im1ymax, im1xmin:im1xmax]
        header = wcsmap.to_header()

        gauss = Gaussian2DKernel(stddev=3)
        im1 = convolve(img, gauss)
        print 'Smoothing finished'

            try:
                fits.writeto('smoothscans/im1_'+region+'.fits', im1, header, clobber=True)

            except IOError:
                os.remove('smoothscans/im1_'+region+'.fits')
                fits.writeto('smoothscans/im1_'+region+'.fits', im1, header, clobber=True)

    if run2 == 1:
        img = img[im1ymin:im1ymax, im1xmin:im1xmax]
        wcsmap = wcsmap[im1ymin:im1ymax, im1xmin:im1xmax]
        bkgd = fits.open('backgrounds/background_im1_'+region+'.fits')[0].data
        im1 = img - bkgd
        header = wcsmap.to_header()

            try:
                fits.writeto('bksubscans/im1_'+region+'.fits', im1, header, clobber=True)

            except IOError:
                os.remove('bksubscans/im1_'+region+'.fits')
                fits.writeto('bksubscans/im1_'+region+'.fits', im1, header, clobber=True)

    print 'im1 saved'

    '''
    #########################################################################
    # Run sextractor
    #########################################################################
    if run1 == 1:
        os.system('sextractor smoothscans/im1_'+region+'.fits -c daofind.sex -CATALOG_NAME sextractor/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME backgrounds/background_im1_'+region+'.fits')


        #os.system('sextractor ../../galexscans/im2_'+region+'_masked.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im2_'+region'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im2_'+region+'.fits')


    # With no background step, subtract background prior to this
    if run2 == 1:
        os.system('sextractor bksubscans/im1_'+region+'.fits -c daofind.sex -CATALOG_NAME sextractor/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0')

    # With weights
    #os.system('sex ../Dunmaps/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_'+region+'fwhm.fits -WEIGHT_IMAGE ../Dunmaps/background/background_im1_'+region+'.fits')

    print 'SExtractor finished'

    #########################################################################
    # Get output from sextractor, convert NUV
    #########################################################################
    im1sex = Table.read('sextractor/sex_im1_'+region+'.fits', format='fits')

    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data['X_IMAGE']+xfac)
    y_new = (data['Y_IMAGE']+yfac)
    nuv = -2.5*np.log10(data['FLUX_AUTO']) + 20.08

    data['x_new'] = x_new
    data['y_new'] = y_new
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

    ascii.write(alldata, 'starcat/starcat_'+currregion+'_01_10_2018.txt', format='ipac')

    print 'Converted to gl, gb, finished'
    '''
