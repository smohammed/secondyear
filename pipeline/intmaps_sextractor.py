from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import os
from astropy.convolution import convolve, Gaussian2DKernel
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# To add scans:
# 1. Add fieldN value and coordinates to skyrange
# 2. Manually check where you want to cut out the images and add to xminmax, yminmax

# Files written:
# 1. Cutouts of main image to remove edges
# 2. SExtractor output, background used
# 3. starcatalog*.fits which adds WCS data to stars
# 4. If yes, sex_galex_matches* includes GALEX matches with 3.5" search radius


#########################################################################
# Select desired field from list
#########################################################################
# All
#skyrange = ['5', '1.4',
skyrange = ['2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

# Regions 7, 8, 9 and 10
#skyrange = ['130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '263.3', '270.5', '288.5', '293.9']

# Test regions
#skyrange = ['4.1', '5.0', '18.5', '19.4', '328.1', '329.0']
#skyrange = ['329.0']
#filt = 'bg_64_filt3'
filt = 'mexhat_4.0_9x9'
code = 'det_thresh2.5_phot_autopar2.5_3.5'

#galex = fits.open('../GALEXAIS.fits')[1].data
#galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')

#########################################################################
# Choose a field as defined above
#########################################################################
for currregion in skyrange:
    #########################################################################
    # Figure out what the cutout range should be
    #########################################################################
    print 'current region = ' + currregion
    region = currregion.replace('.', '')
    img = fits.open('../Dunmaps/count_map_name_'+region+'_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 3532, 51230

    #########################################################################
    # Make cutouts of initial image to help with background correction
    #########################################################################
    img = img[im1ymin:im1ymax, im1xmin:im1xmax]
    gauss = Gaussian2DKernel(stddev=3)
    im1 = convolve(img, gauss)

    print 'Smoothing finished'

    try:
        fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_'+region+'.fits')
    except IOError:
        os.remove('../Dunmaps/im1_'+region+'.fits')
        #os.remove('../Dunmaps/background_im1_'+region+'.fits')

        fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_'+region+'.fits')

    print 'im1 saved'

    #########################################################################
    # Run sextractor
    #########################################################################
    os.system('sex ../Dunmaps/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_'+region+'.fits')

    print 'SExtractor finished'

    #########################################################################
    # Get output from sextractor, convert to gl, gb, NUV
    #########################################################################
    im1sex = Table.read('../Dunmaps/sex_im1_'+region+'.fits', format='fits')
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data['X_IMAGE']+xfac)
    y_new = (data['Y_IMAGE']+yfac)
    nuv = -2.5*np.log10(data['FLUX_AUTO']) + 20.08

    data['x_new'] = x_new
    data['y_new'] = y_new
    data['nuv'] = nuv

    print 'Converted coordinates to gl/gb'

    #########################################################################
    # Combine all tables
    #########################################################################
    #tottable = Table.read('../Dunmaps/sex_im1_'+region+'_edit.fits', format='fits')
    #ascii.write(tottable, '../Dunmaps/sex_total_'+region+'.txt')

    tottable = data

    print 'Combined all data tables'

    #########################################################################
    # Get WCS info
    #########################################################################
    hdulist = fits.open('../Dunmaps/count_map_name_'+region+'_gal_sec_in.fits')

    xpix = tottable['x_new']
    ypix = tottable['y_new']
    w = wcs.WCS(hdulist[0].header)
    pixels = np.array([xpix, ypix]).T
    world = w.wcs_pix2world(pixels, 1)

    glval, gbval = [], []

    for i in range(len(world)):
        glval.append(world[i][0])
        gbval.append(world[i][1])

    if region == '5':
        for i in range(len(glval)):
            if glval[i] > 350:
                glval[i] = glval[i] - 360

    skygal = SkyCoord(glval*u.deg, gbval*u.deg, frame='galactic')
    raval = skygal.icrs.ra.degree
    decval = skygal.icrs.dec.degree

    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcat_'+currregion+'.txt', format='ipac')

    #os.remove('../Dunmaps/im1_'+region+'.fits')

    print 'Added WCS info, finished'

'''
    sex = alldata
    sexgal = SkyCoord(sex['gl']*u.deg, sex['gb']*u.deg, frame='galactic')

    sexind,  galexind,  angsep,  ang3d = search_around_sky(sexgal, galexgal, 3.5*u.arcsec)
    s2 = Table(sex[sexind])
    g2 = Table(galex[galexind])
    comb = hstack([s2, g2])
    comb['angsep'] = angsep
    comb.rename_column('nuv', 'nuv_sex')
    comb.rename_column('nuv_mag', 'nuv_galex')
    comb.rename_column('gl', 'gl_sex')
    comb.rename_column('gb', 'gb_sex')
    comb.rename_column('ra_1', 'ra_sex')
    comb.rename_column('dec_1', 'dec_sex')
    comb.rename_column('ra_2', 'ra_galex')
    comb.rename_column('dec_2', 'dec_galex')
    comb.rename_column('glon', 'gl_galex')
    comb.rename_column('glat', 'gb_galex')
    combcut = np.where(comb['nuv_galex'] == -999.)
    comb.remove_rows(combcut)
    ascii.write(comb, '../Dunmaps/sex_galex_matches_'+region+'_'+filt+code+'.txt', format='basic')

    plt.scatter(comb['nuv_galex'], comb['nuv_sex']-comb['nuv_galex'], alpha=0.1, edgecolor='none')
    plt.axhline(y=0, c='green')
    plt.xlabel('NUV$_{GAIS}$')
    plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
    plt.title('gl ='+str(region)+', '+filt+', '+code)
    plt.xlim((12, 23))
    plt.ylim((-3, 2))
    plt.xlim((12, 23))
    plt.ylim((-3, 2))
    plt.annotate('N = '+str(len(comb)), xy=(13, -2.5))
    plt.savefig('../nuvcomp/03-30-nuvconv_'+region+'_'+filt+code+'.png')
    plt.clf()


#########################################################################
# Now match to other catalogs
#########################################################################
galex = fits.open('../GALEXAIS.fits')[1].data
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')
sex = alldata  # Table.read('../Dunmaps/starcatalog_2057-2102.txt', format='ascii')
sexgal = SkyCoord(sex['gl']*u.deg, sex['gb']*u.deg, frame='galactic')

sexind,  galexind,  angsep,  ang3d = search_around_sky(sexgal, galexgal, 3.5*u.arcsec)
s2 = Table(sex[sexind])
g2 = Table(galex[galexind])

delgl = g2['glon'] - s2['gl']
delgb = g2['glat'] - s2['gb']
dgl = np.mean(g2['glon'] - s2['gl'])
dgb = np.mean(g2['glat'] - s2['gb'])
plt.scatter(delgl*3600, delgb*3600, alpha=0.1)
plt.xlabel('$\Delta$ gl')
plt.ylabel('$\Delta$ gb')
plt.title('len = '+str(len(delgl))+', dgl = '+str(dgl*3600)[:4]+' dgb = '+str(dgb*3600)[:4])
plt.savefig('../Dunmaps/coord_'+region+'.png')
plt.clf()

print 'dgl1 = ', dgl * 3600
print 'dgb1 = ', dgb * 3600

sex['gl'] = sex['gl'] + dgl
sex['gb'] = sex['gb'] + dgb

sexgal = SkyCoord(sex['gl']*u.deg, sex['gb']*u.deg, frame='galactic')
sexind,  galexind,  angsep,  ang3d = search_around_sky(sexgal, galexgal, 3.5*u.arcsec)
s2 = Table(sex[sexind])
g2 = Table(galex[galexind])

delgl = g2['glon'] - s2['gl']
delgb = g2['glat'] - s2['gb']
dgl = np.mean(g2['glon'] - s2['gl'])
dgb = np.mean(g2['glat'] - s2['gb'])

plt.scatter(delgl*3600, delgb*3600, alpha=0.1)
plt.xlabel('$\Delta$ gl')
plt.ylabel('$\Delta$ gb')
plt.title('Fix, len = '+str(len(delgl))+', dgl = '+str(dgl*3600)[:4]+' dgb = '+str(dgb*3600)[:4])
plt.savefig('../Dunmaps/coord_'+region+'_fix.png')
plt.clf()


comb = hstack([s2, g2])
comb['angsep'] = angsep
comb.rename_column('nuv', 'nuv_sex')
comb.rename_column('nuv_mag', 'nuv_galex')
comb.rename_column('gl', 'gl_sex')
comb.rename_column('gb', 'gb_sex')
comb.rename_column('ra_1', 'ra_sex')
comb.rename_column('dec_1', 'dec_sex')
comb.rename_column('ra_2', 'ra_galex')
comb.rename_column('dec_2', 'dec_galex')
comb.rename_column('glon', 'gl_galex')
comb.rename_column('glat', 'gb_galex')

combcut = np.where(comb['nuv_galex'] == -999.)
comb.remove_rows(combcut)
#ascii.write(comb, '../sex_galex_matches_2111-2138.txt', format='basic')
'''
