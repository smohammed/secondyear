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
field1 = '5'
field2 = '8.6-12.2'
field3 = '17.6-19.4'
field4 = '20.3-25.7'
field5 = '205.7-210.2'
field6 = '211.1-213.8'
field7 = '214.7-217.4'
field8 = '218.3-221.0'
field9 = '223.7-226.4'
field10 = '228.2-231.8'
field11 = '1.4'
field12 = '2.3'
field13 = '3.2'
field14 = '4.1'
skyrange = ['5', '8.6-12.2', '17.6-19.4', '20.3-25.7', '205.7-210.2', '211.1-213.8', '214.7-217.4', '218.3-221.0', '223.7-226.4', '228.2-231.8', '1.4']

chosenfield = field11

#########################################################################
# Choose a field as defined above
#########################################################################
region = [x for x in skyrange if chosenfield in x][0].replace('.', '')

#########################################################################
# Figure out what the cutout range should be
#########################################################################
img = fits.open('../Dunmaps/count_map_name'+region+'_gal_sec_in.fits')[0].data

if region == field1.replace('.', ''):  # gl 0.5
    #img = fits.open('../Dunmaps/count_map_name'+region+'_gal_sec_in.fits')[0].data
    #im1xmin, im1xmax, im1ymin, im1ymax = 720, 4350, 1160, 11400
    im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3870, 3572, 51230

if region == field2.replace('.', ''):  # gl 8.6 - 12.2
    im1xmin, im1xmax, im1ymin, im1ymax = 420, 3240, 3500, 13700
    im2xmin, im2xmax, im2ymin, im2ymax = 420, 2060, 1300, 3500
    im3xmin, im3xmax, im3ymin, im3ymax = 2590, 3240, 1400, 3500

if region == field3.replace('.', ''):  # gl 17.6 - 19.4
    im1xmin, im1xmax, im1ymin, im1ymax = 370, 2050, 1350, 13700

if region == field4.replace('.', ''):  # gl 20.3 - 25.7
    im1xmin, im1xmax, im1ymin, im1ymax = 450, 4330, 1400, 13700

if region == field5.replace('.', ''):  # gl 205.7 - 210.2
    im1xmin, im1xmax, im1ymin, im1ymax = 379, 3707, 1370, 13830

if region == field6.replace('.', ''):  # gl 211.1 - 2213.8
    im1xmin, im1xmax, im1ymin, im1ymax = 254, 2538, 1358, 14037

if region == field7.replace('.', ''):  # gl 214.7 - 217.4
    im1xmin, im1xmax, im1ymin, im1ymax = 265, 2569, 1446, 14056

if region == field8.replace('.', ''):  # gl 218.3 - 221.0
    im1xmin, im1xmax, im1ymin, im1ymax = 310, 2650, 1449, 14056

if region == field9.replace('.', ''):  # gl 223.7 - 226.4
    im1xmin, im1xmax, im1ymin, im1ymax = 376, 2056, 1575, 13700
    im2xmin, im2xmax, im2ymin, im2ymax = 2056, 2638, 1575, 12532

if region == field10.replace('.', ''):  # gl 228.2 - 231.8
    im1xmin, im1xmax, im1ymin, im1ymax = 234, 3050, 1554, 14000

if region == field10.replace('.', ''):  # gl 228.2 - 231.8
    im1xmin, im1xmax, im1ymin, im1ymax = 234, 3050, 1554, 14000

if region == field11.replace('.', ''):  # gl 1.4
    im1xmin, im1xmax, im1ymin, im1ymax = 1232, 3960, 3469, 51082


#########################################################################
# Make cutouts of initial image to help with background correction
#########################################################################
im1 = img[im1ymin:im1ymax, im1xmin:im1xmax]
fits.HDUList([fits.PrimaryHDU(im1)]).writeto('../Dunmaps/im1_'+region+'.fits')

if region == field2.replace('.', ''):
    im2 = img[im2ymin:im2ymax, im2xmin:im2xmax]
    fits.HDUList([fits.PrimaryHDU(im2)]).writeto('../Dunmaps/im2_'+region+'.fits')
    im3 = img[im3ymin:im3ymax, im3xmin:im3xmax]
    fits.HDUList([fits.PrimaryHDU(im3)]).writeto('../Dunmaps/im3_'+region+'.fits')

elif region == field9.replace('.', ''):
    im2 = img[im2ymin:im2ymax, im2xmin:im2xmax]
    fits.HDUList([fits.PrimaryHDU(im2)]).writeto('../Dunmaps/im2_'+region+'.fits')

#########################################################################
# Run sextractor
#########################################################################
os.system('sex ../Dunmaps/im1_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im1_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im1_'+region+'.fits')

if region == field2.replace('.', ''):
    os.system('sex ../Dunmaps/im2_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im2_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im2_'+region+'.fits')
    os.system('sex ../Dunmaps/im3_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im3_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im3_'+region+'.fits')

elif region == field9.replace('.', ''):
    os.system('sex ../Dunmaps/im2_'+region+'.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../Dunmaps/sex_im2_'+region+'.fits -CHECKIMAGE_NAME ../Dunmaps/background_im2_'+region+'.fits')

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

if region == field2.replace('.', ''):
    im2sex = fits.open('../Dunmaps/sex_im2_'+region+'.fits')[1].data
    data = im2sex
    xfac = im2xmin
    yfac = im2ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im2_'+region+'_edit.fits')

    im3sex = fits.open('../Dunmaps/sex_im3_'+region+'.fits')[1].data
    data = im3sex
    xfac = im3xmin
    yfac = im3ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im3_'+region+'_edit.fits')

elif region == field9.replace('.', ''):
    im2sex = fits.open('../Dunmaps/sex_im2_'+region+'.fits')[1].data
    data = im2sex
    xfac = im2xmin
    yfac = im2ymin
    x_new = (data.X_IMAGE+xfac)
    y_new = (data.Y_IMAGE+yfac)
    nuv = -2.5*np.log10(data.FLUX_AUTO) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='x_new', format='1E', array=x_new), fits.Column(name='y_new', format='1E', array=y_new), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../Dunmaps/sex_im2_'+region+'_edit.fits')

print 'Converted coordinates to gl/gb'

#########################################################################
# Combine all tables
#########################################################################
if region == field2.replace('.', ''):
    table1 = Table.read('../Dunmaps/sex_im1_'+region+'_edit.fits', format='fits')
    table2 = Table.read('../Dunmaps/sex_im2_'+region+'_edit.fits', format='fits')
    table3 = Table.read('../Dunmaps/sex_im3_'+region+'_edit.fits', format='fits')
    tottable = vstack([table1, table2, table3])
    ascii.write(tottable, '../Dunmaps/sex_total_'+region+'.txt')

elif region == field9.replace('.', ''):
    table1 = Table.read('../Dunmaps/sex_im1_'+region+'_edit.fits', format='fits')
    table2 = Table.read('../Dunmaps/sex_im2_'+region+'_edit.fits', format='fits')
    tottable = vstack([table1, table2])
    ascii.write(tottable, '../Dunmaps/sex_total_'+region+'.txt')

else:
    tottable = Table.read('../Dunmaps/sex_im1_'+region+'_edit.fits', format='fits')
    ascii.write(tottable, '../Dunmaps/sex_total_'+region+'.txt')

print 'Combined all data tables'

#########################################################################
# Get WCS info
#########################################################################
hdulist = fits.open('../Dunmaps/count_map_name'+region+'_gal_sec_in.fits')

xpix = tottable['x_new']
ypix = tottable['y_new']
w = wcs.WCS(hdulist[0].header)
pixels = np.array([xpix, ypix]).T
world = w.wcs_pix2world(pixels, 1)

glval, gbval = [], []

for i in range(len(world)):
    glval.append(world[i][0])
    gbval.append(world[i][1])

if region == field1.replace('.', ''):
    for i in range(len(glval)):
        if glval[i] > 350:
            glval[i] = glval[i] - 360

skygal = SkyCoord(glval*u.deg, gbval*u.deg, frame='galactic')
raval = skygal.icrs.ra.degree
decval = skygal.icrs.dec.degree

coord = Table([glval, gbval, raval, decval], names=('gl', 'gb', 'ra', 'dec'))
alldata = hstack([tottable, coord])
ascii.write(alldata, '../Dunmaps/starcatalog_'+region+'.txt', format='ipac')

print 'Added WCS info, finished'

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
