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
# All complete scans
#scans = ['0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104', '0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194', '0203', '0212', '0221', '0230', '0239', '0248', '0257', '0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356', '0392', '0428', '0437', '0446', '0455', '0464', '0473', '0482', '0491', '0500', '0671', '0689', '0716', '0743', '0752', '0761', '0770', '0779', '0788', '0797', '0806', '0815', '0824', '0833', '0878', '0887', '0896', '0905', '0914', '0923', '0932', '0941', '0950', '0959', '0968', '0977', '0986', '0995', '1004', '1013', '1022', '1031', '1040', '1049', '1058', '1067']

scans = ['216-226']

# Incomplete scans
#incscans = ['9.5', '14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']

# Must redo these for new scans, when available
#incscandict = dict({'9.5': [1214, 3950, 12318, 52037], '14.9': [1214, 3950, 3532, 29730], '79.7': [1214, 3950, 26493, 51563], '90.5': [1214, 3950, 13600, 51230], '91.4': [1214, 3950, 18800, 51600], '103.1': [1214, 3950, 9300, 51600], '104.0': [1214, 3950, 26400, 51600], '122.9': [1214, 3950, 35000, 51600], '127.4': [1214, 3950, 21800, 51300], '223.7': [1214, 3950, 3300, 47000], '273.2': [1214, 3950, 13100, 52300], '283.1': [1214, 3950, 5700, 51600], '289.4': [1214, 3950, 37000, 51500], '306.5': [1214, 3950, 19000, 51500], '309.2': [1214, 3950, 33400, 52400], '324.5': [1214, 3950, 3200, 45000], '329.9': [1214, 3950, 27600, 44500], '338.0': [1214, 3950, 36000, 53000], '339.8': [1605, 3490, 7200, 53000], '342.5': [1214, 3950, 3000, 42900], '343.4': [1214, 3950, 12200, 52600], '345.2': [1214, 3950, 14600, 40800], '348.8': [1214, 3950, 30400, 52600], '349.7': [1214, 3950, 12100, 52600], '350.6': [1214, 3950, 13200, 52600], '351.5': [1214, 3950, 13380, 52600], '352.4': [1214, 3950, 15500, 52700], '353.3': [1214, 3950, 16400, 52700], '354.2': [1214, 3950, 17000, 52700], '355.1': [1214, 3950, 17700, 52700], '356.0': [1214, 3950, 22700, 52700], '357.8': [1214, 3950, 23600, 52700]})


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