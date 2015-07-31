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


img = fits.open('../newfield/count_05-68_gPr_cata_10_corr.fits')[0].data

im1xmin,im1xmax,im1ymin,im1ymax = 2090,4680,1810,4470
im2xmin,im2xmax,im2ymin,im2ymax = 4680,7230,3350,6100
im3xmin,im3xmax,im3ymin,im3ymax = 7230,10100,5010,7570

##################################################
# Make new fields
##################################################
im1 = img[im1ymin:im1ymax,im1xmin:im1xmax]
fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../newfield/im1.fits')

im2 = img[im2ymin:im2ymax,im2xmin:im2xmax]
fits.HDUList([fits.PrimaryHDU(im2)]).writeto('../newfield/im2.fits')

im3 = img[im3ymin:im3ymax,im3xmin:im3xmax]
fits.HDUList([fits.PrimaryHDU(im3)]).writeto('../newfield/im3.fits')

##################################################
# Run sextractor
##################################################
os.system('sex ../newfield/im1.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im1.fits -CHECKIMAGE_NAME ../newfield/background_im1.fits')
os.system('sex ../newfield/im2.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im2.fits -CHECKIMAGE_NAME ../newfield/background_im2.fits')
os.system('sex ../newfield/im3.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../newfield/sex_im3.fits -CHECKIMAGE_NAME ../newfield/background_im3.fits')

##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################
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

print 'Converted coordinates to gl/gb' 

##################################################
# Combine all tables
##################################################
table1 = Table.read('../newfield/sex_im1_edit.fits',format='fits')
table2 = Table.read('../newfield/sex_im2_edit.fits',format='fits')
table3 = Table.read('../newfield/sex_im3_edit.fits',format='fits')

tottable = vstack([table1,table2,table3])

ascii.write(tottable, '../newfield/sex_total_05-68.txt')

print 'Combined all data tables'

##################################################
# Get WCS info
##################################################

hdulist = fits.open('../newfield/count_05-68_gPr_cata_10_corr.fits')
xpix = tottable['x_new']
ypix = tottable['y_new']
w = wcs.WCS(hdulist[0].header)
pixels = np.array([xpix,ypix]).T
world = w.wcs_pix2world(pixels,1)

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

coord = Table([glval,gbval],names=('gl','gb'))

alldata = hstack([tottable,coord])

ascii.write(alldata,'../newfield/starcatalog_05-68.txt')

##################################################
# Now match with Tycho 2/2MASS catalogs
##################################################
sex = Table.read('../newfield/starcatalog_05-68.txt', format='ascii')
t2 = Table.read('../tycho2_2mass_matches.txt', format='ascii')

sexgal = SkyCoord(sex['gl']*u.degree, sex['gb']*u.degree, frame='galactic')
t2gal = SkyCoord(t2['gl_t2']*u.degree, t2['gb_t2']*u.degree, frame='galactic')

t2ind, sexind, angsep, dist3d = search_around_sky(t2gal, sexgal, 10.*u.arcsec)

plt.hist(angsep*3600., bins=100), plt.show()

sex = sex[sexind]
t3 = t2[t2ind]

print len(sex)

combtable = hstack([sex, t3])

ascii.write(combtable, '../starcatalog_05-68_2mass_t2.txt')

print 'Matched with T2 and 2MASS, finished'
print 'Total objects matched =', len(combtable)

