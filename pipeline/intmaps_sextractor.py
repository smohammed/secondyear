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

# Files written:
# 1. rand_deadtime.fits images
# 2. sex*, background
# 3. *_edit.fits
# 4. sex_total*

field1 = 0  # 0.5 - 6.8
field2 = 0  # 17.6 - 19.4
field3 = 0  # 20.3 - 25.7
field4 = 0  # 8.6 - 12.2
field5 = 0  # 205.7 - 210.2
field6 = 0  # 211.1 - 213.8
field7 = 0  # 214.7 - 217.4
field8 = 1  # 218.3 - 221.0


# Field nums are new fields created by Dun Wang ordered here as I processed them

if field1 == 1:  # gl 0.5 - 6.8
    img = fits.open('../Dunmaps/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 720, 4350, 1160, 11400

if field2 == 1:  # gl 17.6 - 19.4
    img = fits.open('../Dunmaps/count_map_name176-194_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 370, 2050, 1350, 13700

if field3 == 1:  # gl 20.3 - 25.7
    img = fits.open('../Dunmaps/count_map_name203-257_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 450, 4330, 1400, 13700

if field4 == 1:  # gl 8.6 - 12.2
    img = fits.open('../Dunmaps/count_map_name86-122_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 420, 3240, 3500, 13700
    im2xmin, im2xmax, im2ymin, im2ymax = 420, 2060, 1300, 3500
    im3xmin, im3xmax, im3ymin, im3ymax = 2590, 3240, 1400, 3500

if field5 == 1:  # gl 205.7 - 210.2
    img = fits.open('../Dunmaps/count_map_name2057-2102_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 379, 3707, 1370, 13830

if field6 == 1:  # gl 211.1 - 2213.8
    img = fits.open('../Dunmaps/count_map_name2111-2138_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 254, 2538, 1358, 14037

if field7 == 1:  # gl 214.7 - 217.4
    img = fits.open('../Dunmaps/count_map_name2147-2174_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 265, 2569, 1446, 14056

if field8 == 1:  # gl 218.3 - 221.0
    img = fits.open('../Dunmaps/count_map_name2183-2210_gal_sec_in.fits')[0].data
    im1xmin, im1xmax, im1ymin, im1ymax = 310, 2650, 1449, 14056


##################################################
# Make new fields
##################################################
if field1 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_05-68.fits')

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

if field6 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_2111-2138.fits')

if field7 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_2147-2174.fits')

if field8 == 1:
    im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
    fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_2183-2210.fits')


##################################################
# Run sextractor
##################################################
if field1 == 1:
    os.system('sex ../Dunmaps/im1_05-68.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_05-68.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_05-68.fits')

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

if field6 == 1:
    os.system('sex ../Dunmaps/im1_2111-2138.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_2111-2138.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_2111-2138.fits')

if field7 == 1:
    os.system('sex ../Dunmaps/im1_2147-2174.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_2147-2174.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_2147-2174.fits')

if field8 == 1:
    os.system('sex ../Dunmaps/im1_2183-2210.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_2183-2210.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_2183-2210.fits')


##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################
if field1 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_05-68.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_05-68_edit.fits')

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

if field6 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_2111-2138.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_2111-2138_edit.fits')

if field7 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_2147-2174.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_2147-2174_edit.fits')

if field8 == 1:
    im1sex = fits.open('../Dunmaps/sex_im1_2183-2210.fits')[1].data
    data = im1sex
    xfac = im1xmin
    yfac = im1ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im1_2183-2210_edit.fits')


print 'Converted coordinates to gl/gb'

##################################################
# Combine all tables
##################################################
if field1 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_05-68_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_05-68.txt')

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

if field6 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_2111-2138_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_2111-2138.txt')

if field7 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_2147-2174_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_2147-2174.txt')

if field8 == 1:
    tottable = Table.read('../Dunmaps/sex_im1_2183-2210_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_2183-2210.txt')


print 'Combined all data tables'

##################################################
# Get WCS info
##################################################
if field1 == 1:
    hdulist = fits.open('../Dunmaps/count_map_05-68_gPr_cata_10_corr_gal.fits')

if field2 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name176-194_gal_sec_in.fits')

if field3 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name203-257_gal_sec_in.fits')

if field4 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name86-122_gal_sec_in.fits')

if field5 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name2057-2102_gal_sec_in.fits')

if field6 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name2111-2138_gal_sec_in.fits')

if field7 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name2147-2174_gal_sec_in.fits')

if field8 == 1:
    hdulist = fits.open('../Dunmaps/count_map_name2183-2210_gal_sec_in.fits')


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

if (field2 == 1) or (field3 == 1) or (field4 == 1) or (field5 == 1) or (field6 == 1) or (field7 == 1) or (field8 == 1):
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
    ascii.write(alldata, '../Dunmaps/starcatalog_05-68.txt', format='ipac')

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

if field6 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_2111-2138.txt', format='ipac')

if field7 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_2147-2174.txt', format='ipac')

if field8 == 1:
    coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
    alldata = hstack([tottable, coord])
    ascii.write(alldata, '../Dunmaps/starcatalog_2183-2210.txt', format='ipac')



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
plt.title('len = '+str(len(delgl))+', dgl = '+str(dgl*3600)[:6]+' dgb = '+str(dgb*3600)[:6])
plt.savefig('../Dunmaps/coord_2183-2210.png')
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
plt.title('Fix, len = '+str(len(delgl))+', dgl = '+str(dgl*3600)[:6]+' dgb = '+str(dgb*3600)[:6])
plt.savefig('../Dunmaps/coord_2183-2210_fix.png')
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

ascii.write(comb, '../sex_galex_matches_2183-2210.txt', format='basic')



'''
a = Table.read('../sex_galex_matches_176-194.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (17.6-19.4)')
plt.savefig('02-10-nuvcomp_sextractor_gais_176-194.png')
plt.clf()

a = Table.read('../sex_galex_matches_203-257.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (20.3-25.7)')
plt.savefig('02-10-nuvcomp_sextractor_gais_203-257.png')
plt.clf()

a = Table.read('../sex_galex_matches_86-122.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (8.6-12.2)')
plt.savefig('02-10-nuvcomp_sextractor_gais_86-122.png')
plt.clf()

a = Table.read('../sex_galex_matches_2057-2102.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (205.7-210.2)')
plt.savefig('02-10-nuvcomp_sextractor_gais_2057-2102.png')
plt.clf()

a = Table.read('../sex_galex_matches_2111-2138.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (211.1-213.8)')
plt.savefig('02-10-nuvcomp_sextractor_gais_2111-2138.png')
plt.clf()

a = Table.read('../sex_galex_matches_2147-2174.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (214.7-217.4)')
plt.savefig('02-10-nuvcomp_sextractor_gais_2147-2174.png')
plt.clf()

a = Table.read('../sex_galex_matches_2183-2210.txt', format='ascii')
plt.scatter(a['nuv_galex'],a['nuv_sex']-a['nuv_galex'],alpha=0.1)
plt.axhline(y=0,c='black')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{GAIS}$ - NUV$_{GAIS}$')
plt.xlim((11,22))
plt.ylim((-2,2))
plt.title('gl = (218.3-221.0)')
plt.savefig('02-10-nuvcomp_sextractor_gais_2183-2210.png')
plt.clf()
'''