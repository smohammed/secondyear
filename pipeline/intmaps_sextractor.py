from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import os
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# To add scans:
# 1. Add fieldN value and coordinates to skyrange
# 2. Manually check where you want to cut out the images and add to xminmax, yminmax

# Files written:
# 1. Cutouts of main image to remove edges
# 2. SExtractor output, background used
# 3. *_edit.fits file that converts pixels back to original image frame and NUV calculation
# 4. sex_total*.fits which combines all tables if there are multiple cutouts
# 5. starcatalog*.fits which adds WCS data to stars
# 6. If yes, sex_galex_matches* includes GALEX matches with 3.5" search radius


#########################################################################
# Select desired field from list
#########################################################################
# Set 1 and 2
#skyrange = ['10.4', '11.3', '12.2', '14.0', '14.9', '1.4', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '2.3', '24.8', '25.7', '28.4', '302.9', '303.8', '304.7', '306.5', '309.2', '316.4', '320.0', '320.9', '32.0', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '328.1', '329.0', '329.9', '32.9', '3.2', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '33.8', '339.8', '341.6', '342.5', '343.4', '345.2', '34.7', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '35.6', '357.8', '358.7', '359.6', '39.2', '4.1', '5.0', '5.9', '5', '6.8', '9.5']

# Set 3
#skyrange = ['290.3', '291.2', '292.1', '293.0', '29.3', '295.7', '297.5', '298.4', '301.1', '302.0', '30.2', '305.6', '308.3', '310.1', '31.1', '315.5', '317.3', '318.2', '319.1', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '8.6', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5']

# Set 4
skyrange = ['102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '271.4', '272.3', '273.2', '274.1', '278.6', '279.5', '281.3', '283.1', '284.0', '289.4', '327.2']

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
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_'+region+'.fits')

    #########################################################################
    # Run sextractor
    #########################################################################
    os.system('sex ../Dunmaps/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_'+region+'.fits')

    print 'SExtractor ran'

    #########################################################################
    # Get output from sextractor, convert to gl, gb, NUV
    #########################################################################
    im1sex = fits.open('../Dunmaps/sex_im1_'+region+'.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_'+region+'_edit.fits')

    print 'Converted coordinates to gl/gb'

    #########################################################################
    # Combine all tables
    #########################################################################
    tottable = Table.read('../Dunmaps/sex_im1_'+region+'_edit.fits', format='fits')
    #ascii.write(tottable, '../Dunmaps/sex_total_'+region+'.txt')

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

    print 'Added WCS info, finished'

'''
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