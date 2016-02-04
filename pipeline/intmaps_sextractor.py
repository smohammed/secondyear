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

field1 = 0  # 0-7
field2 = 0  # 17.6-19.4
field3 = 0  # 20.3-25.7
field4 = 0  # 8.6-12.2
field5 = 1  # 205.7-210.2

if field1 == 1: # gl 0-7
    img = fits.open('../newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 720, 4350, 1160, 11400 # gl 0-7

# Field#s are new fields created by Dun Wang ordered here as I processed them

if field2 == 1: # gl 17.6 - 19.4
    img = fits.open('../Dunmaps/count_map_name176-194_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 370, 2050, 1350, 13700

if field3 == 1: # gl 20.3 - 25.7
    img = fits.open('../Dunmaps/count_map_name203-257_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 450, 4330, 1400, 13700

if field4 == 1: # gl 8.6 - 12.2
    img = fits.open('../Dunmaps/count_map_name86-122_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 420, 3240, 3500, 13700
    im2xmin, im2xmax, im2ymin, im2ymax = 420, 2060, 1300, 3500
    im3xmin, im3xmax, im3ymin, im3ymax = 2590, 3240, 1400, 3500

if field5 == 1: # gl 205.7 - 210.2
    img = fits.open('../Dunmaps/count_map_name2057-2102_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 379,3707,1370,13830


##################################################
# Make new fields
##################################################
if field1 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../newfield/im1_gal_tests_2.fits')

if field2 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_176-194.fits')

if field3 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_203-257.fits')

if field4 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_86-122.fits')

    im2 = img[im2ymin:im2ymax, im2xmin:im2xmax]
    fits.HDUList([fits.PrimaryHDU(im2)]).writeto('../Dunmaps/im2_86-122.fits')

    im3 = img[im3ymin:im3ymax, im3xmin:im3xmax]
    fits.HDUList([fits.PrimaryHDU(im3)]).writeto('../Dunmaps/im3_86-122.fits')

if field5 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_2057-2102.fits')


##################################################
# Run sextractor
##################################################
if field1 == 1:
    os.system('sex ../newfield/im1_gal_tests_2.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im1_gal_tests_2.fits -CHECKIMAGE_NAME ../newfield/background_im1_gal_tests_2.fits')

if field2 == 1:
    os.system('sex ../Dunmaps/im1_176-194.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_176-194.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_176-194.fits')

if field3 == 1:
    os.system('sex ../Dunmaps/im1_203-257.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_203-257.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_203-257.fits')

if field4 == 1:
    os.system('sex ../Dunmaps/im1_86-122.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_86-122.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_86-122.fits')
    os.system('sex ../Dunmaps/im2_86-122.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im2_86-122.fits -CHECKIMAGE_NAME ../Dunmaps/background_im2_86-122.fits')
    os.system('sex ../Dunmaps/im3_86-122.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im3_86-122.fits -CHECKIMAGE_NAME ../Dunmaps/background_im3_86-122.fits')

if field5 == 1:
    os.system('sex ../Dunmaps/im1_2057-2102.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_2057-2102.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_2057-2102.fits')


##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################
if field1 == 1:
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

if field2 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_176-194.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_176-194_edit.fits')

if field3 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_203-257.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_203-257_edit.fits')

if field4 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_86-122.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_86-122_edit.fits')

    im2sex = fits.open('../Dunmaps/sex_im2_86-122.fits')[1].data
    data = im2sex
    xfac = im2xmin
    yfac = im2ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im2_86-122_edit.fits')

    im3sex = fits.open('../Dunmaps/sex_im3_86-122.fits')[1].data
    data = im3sex
    xfac = im3xmin
    yfac = im3ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im3_86-122_edit.fits')

if field5 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_2057-2102.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_2057-2102_edit.fits')

print 'Converted coordinates to gl/gb'

##################################################
# Combine all tables
##################################################
if field1 == 1:
    tottable = Table.read('../newfield/sex_im1_gal_tests_2_edit.fits', format='fits')
    ascii.write(tottable, '../newfield/sex_total_05-68_gal_tests_2.txt')

if field2 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_176-194_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_176-194.txt')

if field3 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_203-257_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_203-257.txt')

if field4 == 1:
    table1 = Table.read('../Dunmaps/sex_im1_86-122_edit.fits', format='fits')
    table2 = Table.read('../Dunmaps/sex_im2_86-122_edit.fits', format='fits')
    table3 = Table.read('../Dunmaps/sex_im3_86-122_edit.fits', format='fits')
    tottable = vstack([table1, table2, table3])
    ascii.write(tottable, '../Dunmaps/sex_total_86-122.txt')

if field5 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_2057-2102_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_2057-2102.txt')

print 'Combined all data tables'

##################################################
# Get WCS info
##################################################
if field1 == 1:
    hdulist = fits.open('../newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')

if field2 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name176-194_gal_sec_in.fits')

if field3 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name203-257_gal_sec_in.fits')

if field4 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name86-122_gal_sec_in.fits')

if field5 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name2057-2102_gal_sec_in.fits')


xpix = tottable['x_new']
ypix = tottable['y_new']
w = wcs.WCS(hdulist[0].header)
pixels = np.array([xpix, ypix]).T
world = w.wcs_pix2world(pixels, 1)

if field1 == 1:
    glval = []
    gbval = []

    for i in range(len(world)):
        glval.append(world[i][0])
        gbval.append(world[i][1])

    for i in range(len(glval)):
        if glval[i] > 350:
            glval[i] = glval[i] - 360

    skygal = SkyCoord(glval*u.deg, gbval*u.deg, frame='galactic')
    raval = skygal.icrs.ra.degree
    decval = skygal.icrs.dec.degree

if (field2 == 1) or (field3 == 1) or (field4 == 1) or (field5 == 1):
    glval = []
    gbval = []

    for i in range(len(world)):
        glval.append(world[i][0])
        gbval.append(world[i][1])

    skygal = SkyCoord(glval*u.deg, gbval*u.deg, frame='galactic')
    raval = skygal.icrs.ra.degree
    decval = skygal.icrs.dec.degree


if field1 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    '''
    alldata.rename_column('cntr_01','cntr')
    alldata.rename_column('number_01','number')
    alldata.rename_column('x_image_01','x_image')
    alldata.rename_column('y_image_01','y_image')
    alldata.rename_column('flux_auto_01','flux_auto')
    alldata.rename_column('fluxerr_auto_01','fluxerr_auto')
    alldata.rename_column('x_new_01','x_new')
    alldata.rename_column('y_new_01','y_new')
    alldata.rename_column('nuv_01','nuv')
    alldata.rename_column('gl_01','gl_sex')
    alldata.rename_column('gb_01','gb_sex')
    alldata.rename_column('ra_01','ra_sex')
    alldata.rename_column('dec_01','dec_sex')
    alldata.rename_column('ra','ra_2mass')
    alldata.rename_column('dec','dec_2mass')
    alldata.rename_column('j_m','j')
    alldata.rename_column('h_m','h')
    alldata.rename_column('k_m','k')
    '''
    ascii.write(alldata, '../newfield/starcatalog_05-68_gal_tests_2.txt', format='ipac')

if field2 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_176-194.txt', format='ipac')

if field3 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_203-257.txt', format='ipac')

if field4 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_86-122.txt', format='ipac')

if field5 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_2057-2102.txt', format='ipac')


galex = fits.open('../GALEXAIS.fits')[1].data
galexgal = SkyCoord(galex['ra']*u.deg, galex['dec']*u.deg, frame='icrs')
a = Table.read('../Dunmaps/starcatalog_2057-2102.txt', format='ascii')
agal = SkyCoord(a['ra']*u.deg, a['dec']*u.deg, frame='icrs')

aind,  galexind,  angsep,  ang3d = search_around_sky(agal, galexgal, 3.5*u.arcsec)
a2 = Table(a[aind])
g2 = Table(galex[galexind])
comb = hstack([a2, g2])
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
ascii.write(comb, '../sex_galex_matches_2057-2102.txt', format='basic')
