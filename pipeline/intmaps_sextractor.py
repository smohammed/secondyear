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
#scans = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '10.4', '11.3', '12.2', '14.0', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '102.2', '104.9', '105.8', '106.7', '107.6', '110.3', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '124.7', '125.6', '126.5', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '284.0', '285.8', '286.7', '288.5', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '308.3', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '325.4', '326.3', '327.2', '328.1', '329.0', '331.7', '332.6', '333.5', '334.4', '335.3', '338.9', '341.6', '358.7', '359.6']

# fec scans
#scans = ['0014', '0032', '0059', '0203', '0239', '0356', '0392', '0743', '1103']
#scans = ['0014']#, '0032','0059', '0203', '0239']  # These scans supposively has the half pixel fix, as of 04/17

#scans = ['0005', '0014', '0023','0032','0041','0050','0059','0068','0086','0095','0104']

#scans = ['0014', '0086']

#scans = ['0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104',

scans = ['0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194', '0203', '0212', '0221', '0230', '0239', '0248', '0257', '0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356', '0392', '0428', '0437', '0446', '0455', '0464', '0473', '0482', '0491', '0500', '0671', '0689', '0716', '0743', '0752', '0761', '0770', '0779', '0788', '0797', '0806', '0815', '0824', '0833', '0878', '0887', '0896', '0905', '0914', '0923', '0932', '0941', '0950', '0959', '0968', '0977', '0986', '0995', '1004', '1013', '1022', '1031', '1040', '1049', '1058', '1067']


#scans = ['0113', '0113', '0113', '0122', '0122', '0122', '0140', '0140', '0140', '0149', '0149', '0149', '0158', '0158', '0158', '0167', '0167', '0167', '0176', '0176', '0176', '0185', '0185', '0185', '0194', '0194', '0194', '0203', '0203', '0203', '0212', '0212', '0212', '0221', '0221', '0221', '0230', '0230', '0230', '0239', '0239', '0239', '0248', '0248', '0248', '0257', '0257', '0257', '0284', '0284', '0284', '0293', '0293', '0293', '0302', '0302', '0302', '0311', '0311', '0311', '0320', '0320', '0320', '0329', '0329', '0329', '0338', '0338', '0338', '0347', '0347', '0347', '0356', '0356', '0356', '0392', '0392', '0392', '0428', '0428', '0428', '0437', '0437', '0437', '0446', '0446', '0446', '0455', '0455', '0455', '0464', '0464', '0464', '0473', '0473', '0473', '0482', '0482', '0482', '0491', '0491', '0491', '0500', '0500', '0500']

# Incomplete scans
incscans = ['9.5', '14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']

# Must redo these for new scans, when available
incscandict = dict({'9.5': [1214, 3950, 12318, 52037], '14.9': [1214, 3950, 3532, 29730], '79.7': [1214, 3950, 26493, 51563], '90.5': [1214, 3950, 13600, 51230], '91.4': [1214, 3950, 18800, 51600], '103.1': [1214, 3950, 9300, 51600], '104.0': [1214, 3950, 26400, 51600], '122.9': [1214, 3950, 35000, 51600], '127.4': [1214, 3950, 21800, 51300], '223.7': [1214, 3950, 3300, 47000], '273.2': [1214, 3950, 13100, 52300], '283.1': [1214, 3950, 5700, 51600], '289.4': [1214, 3950, 37000, 51500], '306.5': [1214, 3950, 19000, 51500], '309.2': [1214, 3950, 33400, 52400], '324.5': [1214, 3950, 3200, 45000], '329.9': [1214, 3950, 27600, 44500], '338.0': [1214, 3950, 36000, 53000], '339.8': [1605, 3490, 7200, 53000], '342.5': [1214, 3950, 3000, 42900], '343.4': [1214, 3950, 12200, 52600], '345.2': [1214, 3950, 14600, 40800], '348.8': [1214, 3950, 30400, 52600], '349.7': [1214, 3950, 12100, 52600], '350.6': [1214, 3950, 13200, 52600], '351.5': [1214, 3950, 13380, 52600], '352.4': [1214, 3950, 15500, 52700], '353.3': [1214, 3950, 16400, 52700], '354.2': [1214, 3950, 17000, 52700], '355.1': [1214, 3950, 17700, 52700], '356.0': [1214, 3950, 22700, 52700], '357.8': [1214, 3950, 23600, 52700]})


# Custom parameters that I use
#filt = 'mexhat_4.0_9x9'
filt = 'gauss_3.0_7x7.conv'
code = 'det_thresh4_phot_autopar2.5_3.5'
# Det Thresh = 4, Analysis Thresh = 3.5

def createCircularMask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask


#########################################################################
# Decide to run on full or partial scans
#########################################################################
run1 = 1
run2 = 0
run3 = 0

full = 1
partial = 0

if full == 1:
    skyrange = scans
if partial == 1:
    skyrange = incscans

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
    if full == 1:
        #im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 3532, 51230 # Old range
        im1xmin, im1xmax, im1ymin, im1ymax = 379, 2515, 2000, 38000 # Out of 
    if partial == 1:
        im1xmin, im1xmax, im1ymin, im1ymax = incscandict[currregion][0], incscandict[currregion][1], incscandict[currregion][2], incscandict[currregion][3]
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

    if run2 == 1:
        im2 = img.copy()
        catalog = Table.read('../../galexscans/starcat_'+currregion+'_11-10.txt', format='ascii')
        c2 = catalog[np.where(catalog['nuv'] < 13)]
        h, w = img.shape[:2]

        for i in range(len(c2)):
            mask = createCircularMask(h, w, center=[c2['X_IMAGE'][i], c2['Y_IMAGE'][i]], radius=50)
            im2[mask] = 0

        fits.writeto('../../galexscans/im2_'+region+'_masked.fits', im2, hdu.header, clobber=True)
        
        
    if run3 == 1:
        img = img[im1ymin:im1ymax, im1xmin:im1xmax]
        wcsmap = wcsmap[im1ymin:im1ymax, im1xmin:im1xmax]
        bkgd = fits.open('../../galexscans/background_im2_'+region+'.fits')[0].data
        im1 = img - bkgd
        header = wcsmap.to_header()

    if (run1 or run3):
        try:
            fits.writeto('../../galexscans/im1_'+region+'.fits', im1, header, clobber=True)

        except IOError:
            os.remove('../../galexscans/im1_'+region+'.fits')
            fits.writeto('../../galexscans/im1_'+region+'.fits', im1, header, clobber=True)

    print 'im1 saved'

    #########################################################################
    # Run sextractor
    #########################################################################
    if run1 == 1:
        os.system('sextractor ../../galexscans/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im1_'+region+'.fits')


    # Now on masked image, masking NUV < 13
    if run2 == 1:
        os.system('sextractor ../../galexscans/im2_'+region+'_masked.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im2_'+region+'.fits')

        #os.system('sextractor ../../galexscans/im2_'+region+'_masked.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im2_'+region'.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME ../../galexscans/background_im2_'+region+'.fits')


    # With no background step, subtract background prior to this
    if run3 == 1:
        os.system('sextractor ../../galexscans/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../../galexscans/sex_im1_'+region+'.fits -BACK_TYPE MANUAL -BACK_VALUE 0.0')

    # With weights
    #os.system('sex ../Dunmaps/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_'+region+'fwhm.fits -WEIGHT_IMAGE ../Dunmaps/background/background_im1_'+region+'.fits')

    print 'SExtractor finished'

    #########################################################################
    # Get output from sextractor, convert NUV
    #########################################################################
    im1sex = Table.read('../../galexscans/sex_im1_'+region+'.fits', format='fits')

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

    ascii.write(alldata, '../../galexscans/starcat_'+currregion+'_11-10.txt', format='ipac')

    print 'Converted to gl, gb, finished'
