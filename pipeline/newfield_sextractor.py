from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import sys
import os
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# Files written:
# 1. rand_deadtime.fits images
# 2. sex*, background
# 3. *_edit.fits
# 4. sex_total*
# 5. t2_2mass.fits match

olddata = 0
newdata = 1

if olddata == 1:
    img = fits.open('../newfield/count_05-68_gPr_cata_10_corr.fits')[0].data
    im1xmin,im1xmax,im1ymin,im1ymax = 2090,4680,1810,4470
    im2xmin,im2xmax,im2ymin,im2ymax = 4680,7230,3350,6100
    im3xmin,im3xmax,im3ymin,im3ymax = 7230,10100,5010,7570

if newdata == 1:
    img = fits.open('../newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
    im1xmin,im1xmax,im1ymin,im1ymax = 720,4350,1160,11400

##################################################
# Make new fields
##################################################
if olddata == 1:
    im1 = img[im1ymin:im1ymax,im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../newfield/im1.fits')
    
    im2 = img[im2ymin:im2ymax,im2xmin:im2xmax]
    fits.HDUList([fits.PrimaryHDU(im2)]).writeto('../newfield/im2.fits')
    
    im3 = img[im3ymin:im3ymax,im3xmin:im3xmax]
    fits.HDUList([fits.PrimaryHDU(im3)]).writeto('../newfield/im3.fits')

if newdata == 1:
    im1 = img[im1ymin:im1ymax,im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../newfield/im1_gal_tests_2.fits')


##################################################
# Run sextractor
##################################################
if olddata == 1:
    os.system('sex ../newfield/im1.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im1.fits -CHECKIMAGE_NAME ../newfield/background_im1.fits')
    os.system('sex ../newfield/im2.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im2.fits -CHECKIMAGE_NAME ../newfield/background_im2.fits')
    os.system('sex ../newfield/im3.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im3.fits -CHECKIMAGE_NAME ../newfield/background_im3.fits')

if newdata == 1:
    os.system('sex ../newfield/im1_gal_tests_2.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im1_gal_tests_2.fits -CHECKIMAGE_NAME ../newfield/background_im1_gal_tests_2.fits')


##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################
if olddata == 1:
    im1sex = fits.open('../newfield/sex_im1.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac) 
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../newfield/sex_im1_edit.fits')
    
    im2sex = fits.open('../newfield/sex_im2.fits')[1].data
    data = im2sex
    xfac = im2xmin
    yfac = im2ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac) 
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../newfield/sex_im2_edit.fits')
    
    im3sex = fits.open('../newfield/sex_im3.fits')[1].data
    data = im3sex
    xfac = im3xmin
    yfac = im3ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac) 
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../newfield/sex_im3_edit.fits')

if newdata == 1:
    im1sex = fits.open('../newfield/sex_im1_gal_tests_2.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac) 
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../newfield/sex_im1_gal_tests_2_edit.fits')


print 'Converted coordinates to gl/gb' 

##################################################
# Combine all tables
##################################################
if olddata == 1:
    table1 = Table.read('../newfield/sex_im1_edit.fits',format='fits')
    table2 = Table.read('../newfield/sex_im2_edit.fits',format='fits')
    table3 = Table.read('../newfield/sex_im3_edit.fits',format='fits')
    tottable = vstack([table1,table2,table3])
    ascii.write(tottable, '../newfield/sex_total_05-68.txt')

if newdata == 1:
    tottable = Table.read('../newfield/sex_im1_gal_tests_2_edit.fits',format='fits')
    ascii.write(tottable, '../newfield/sex_total_05-68_gal_tests_2.txt')    

print 'Combined all data tables'

##################################################
# Get WCS info
##################################################
if olddata == 1:
    hdulist = fits.open('../newfield/count_05-68_gPr_cata_10_corr.fits')

if newdata == 1:
    hdulist = fits.open('../newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')

xpix = tottable['x_new']
ypix = tottable['y_new']
w = wcs.WCS(hdulist[0].header)
pixels = np.array([xpix,ypix]).T
world = w.wcs_pix2world(pixels,1)

if olddata == 1:
    raval = []
    decval = []
    
    for i in range(len(world)):
        raval.append(world[i][0])
        decval.append(world[i][1])
    
    skygal = SkyCoord(raval*u.deg,decval*u.deg,frame='icrs')
    glval = skygal.galactic.l.degree
    gbval = skygal.galactic.b.degree
    
    for i in range(len(glval)):
        if glval[i] > 350:
            glval[i] = glval[i] - 360

if newdata == 1:
    glval = []
    gbval = []
    
    for i in range(len(world)):
        glval.append(world[i][0])
        gbval.append(world[i][1])

    for i in range(len(glval)):
        if glval[i] > 350:
            glval[i] = glval[i] - 360

    skygal = SkyCoord(glval*u.deg,gbval*u.deg,frame='galactic')
    raval = skygal.icrs.ra.degree
    decval = skygal.icrs.dec.degree


if olddata == 1:
    coord = Table([glval,gbval],names=('gl','gb'))
    alldata = hstack([tottable,coord])
    ascii.write(alldata,'../newfield/starcatalog_05-68.txt')

if newdata == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../newfield/starcatalog_05-68_gal_tests_2.txt',format='ipac')



# Then manually match with 2MASS
# Use J < 13.5 mag and 3 arcsec search

##################################################
# Now match with Tycho 2/2MASS catalogs
##################################################

sex = Table.read('newfield_ipac_gal_2mass_matches_3arcsec_J_lt13.5.txt', format='ascii')
tycho = Table.read('tycho2.fits', format='fits')

tcut = np.where((tycho['Glon'] > 0.) & (tycho['Glon'] < 8) & (tycho['Glat'] > -10) & (tycho['Glat'] < 10))

t2 = tycho[tcut]

sexgal = SkyCoord(sex['ra_2mass']*u.degree, sex['dec_2mass']*u.degree, frame='icrs')
t2gal = SkyCoord(t2['RAJ2000']*u.degree, t2['DEJ2000']*u.degree, frame='icrs')

t2ind, sexind, angsep, dist3d = search_around_sky(t2gal, sexgal, 1.*u.arcsec)

plt.hist(angsep*3600., bins=100), plt.show()

s2 = sex[sexind]
t3 = t2[t2ind]

print len(s2)

combtable = hstack([s2, t3])

ascii.write(combtable, '../newfield_gal_2mass_t2_jlim_13.5_3arcsec_tests.txt')

print 'Matched with T2 and 2MASS, finished'
print 'Total objects matched =', len(combtable)

