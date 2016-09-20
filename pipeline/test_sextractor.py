from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import os
from astropy.convolution import convolve, Gaussian2DKernel
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

galex = fits.open('../GALEXAIS.fits')[1].data
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')


region = '3497'

smooth = 0

#########################################################################
# Figure out what the cutout range should be
#########################################################################
img = fits.open('../Dunmaps/countmaps/count_map_name_'+region+'_gal_sec_in.fits')[0].data
#im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 3532, 51230
im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 12100, 52600

#########################################################################
# Make cutouts of initial image to help with background correction
#########################################################################
img = img[im1ymin:im1ymax, im1xmin:im1xmax]

im1 = img + 0.01

#print 'Loaded image'

# Smooth to get background map
if smooth == 1:
    gauss = Gaussian2DKernel(stddev=3)
    im1 = convolve(img, gauss)
    print 'Smoothing finished'

else:
    bkgd = fits.open('background_im1_'+region+'_tests.fits')[0].data
    im1 = img - bkgd
    expmap = fits.open('../../galexfiles/exp/count_map_name_'+region+'_gal_sec_exp.fits')[0].data
    expmap = expmap[im1ymin:im1ymax, im1xmin:im1xmax]

try:
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('im1_'+region+'_tests.fits')
except IOError:
    os.remove('im1_'+region+'_tests.fits')
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('im1_'+region+'_tests.fits')

if smooth == 0:
    try:
        fits.HDUList([fits.PrimaryHDU(expmap)]).writeto('exp_'+region+'_tests.fits')
    except IOError:
        os.remove('exp_'+region+'_tests.fits')
        fits.HDUList([fits.PrimaryHDU(expmap)]).writeto('exp_'+region+'_tests.fits')

print 'im1 saved'

#########################################################################
# Run sextractor
#########################################################################
if smooth == 1:
    os.system('sex im1_'+region+'_tests.fits -c ~/sextractor/daofind.sex -CATALOG_NAME sex_im1_'+region+'_tests.fits -BACK_TYPE AUTO -CHECKIMAGE_NAME background_im1_'+region+'_tests.fits')

else:
    # No weights
    #os.system('sex im1_'+region+'_tests.fits -c ~/sextractor/daofind.sex -CATALOG_NAME sex_im1_'+region+'_tests.fits -BACK_TYPE AUTO')

    # EXP weight
    os.system('sex im1_'+region+'_tests.fits -c ~/sextractor/daofind.sex -CATALOG_NAME sex_im1_'+region+'_tests.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE exp_'+region+'_tests.fits')

# With weights
#os.system('sex im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME sex_im1_'+region+'fwhm_nuvtest.fits -WEIGHT_IMAGE background_im1_'+region+'_smooth.fits') #../Dunmaps/background/background_im1_'+region+'.fits')


print 'SExtractor finished'

#########################################################################
# Get output from sextractor, convert to gl, gb, NUV
#########################################################################
im1sex = Table.read('sex_im1_'+region+'_tests.fits', format='fits')
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

print 'Converted coordinates to gl/gb'

#########################################################################
# Get WCS info
#########################################################################
hdulist = fits.open('../Dunmaps/countmaps/count_map_name_'+region+'_gal_sec_in.fits')

xpix = data['x_new']
ypix = data['y_new']
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
alldata = hstack([data, coord])

alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5))]
ascii.write(alldata, 'starcat_'+region+'_tests.txt', format='ipac')

os.remove('im1_'+region+'_tests.fits')
print 'Added WCS info, finished'

a = alldata
agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
aind,  gaind,  angsepa,  ang3d = search_around_sky(agal, galexgal, 3.5*u.arcsec)
a2 = Table(a[aind])
g2 = Table(galex[gaind])
comb = hstack([a2, g2])
comb['angsep'] = angsepa
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

plt.scatter(comb['nuv_galex'], comb['nuv_sex']-comb['nuv_galex'], edgecolor='none', s=5, c=comb['gb_sex'], vmin=-10, vmax=10)
plt.colorbar()
plt.axhline(y=0, c='red')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.title('gl = '+region)
plt.xlim((13, 25))
plt.ylim((-1, 1.5))
plt.annotate('N = '+str(len(comb)), xy=(22, 1))
plt.show()
#ascii.write(comb, 'sex_galex_matches_'+region+'_nuvtest.txt', format='basic')

print 'GALEX matched'
