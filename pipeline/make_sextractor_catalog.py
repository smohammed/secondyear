from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
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


glstart = 20

img = fits.open('../combmaps12000/comb12000_gl'+str(glstart)+'to'+str(glstart+20)+'.fits')[0].data

if glstart == 20:
    im1xmin,im1xmax,im1ymin,im1ymax = 0,3770,0,11999
    im2xmin,im2xmax,im2ymin,im2ymax = 4750,5370,0,11999
    im3xmin,im3xmax,im3ymin,im3ymax = 5370,6340,7900,11999
    im4xmin,im4xmax,im4ymin,im4ymax = 6340,9700,0,11999
    im5xmin,im5xmax,im5ymin,im5ymax = 11230,11870,0,11999

if glstart == 40:
    im1xmin,im1xmax,im1ymin,im1ymax = 1380,6350,0,11999

if glstart == 60:
    im1xmin,im1xmax,im1ymin,im1ymax = 3950,4600,0,11999
    im2xmin,im2xmax,im2ymin,im2ymax = 5020,5700,0,11999
    im3xmin,im3xmax,im3ymin,im3ymax = 6640,7320,0,11999
    im4xmin,im4xmax,im4ymin,im4ymax = 8260,11999,11720,11999
    im5xmin,im5xmax,im5ymin,im5ymax = 8260,11620,0,11720

if glstart == 80:
    im1xmin,im1xmax,im1ymin,im1ymax = 30,2330,50,11999
    im2xmin,im2xmax,im2ymin,im2ymax = 4400,7060,3800,11999
    im3xmin,im3xmax,im3ymin,im3ymax = 4400,6070,0,3800
    im4xmin,im4xmax,im4ymin,im4ymax = 6070,6620,2600,3800
    im5xmin,im5xmax,im5ymin,im5ymax = 7060,11999,0,11999

if glstart == 100: # CHECK IM6 HERE BECAUSE IT DOESN'T BELONG HERE
    im1xmin,im1xmax,im1ymin,im1ymax = 0,1670,0,5900
    im2xmin,im2xmax,im2ymin,im2ymax = 1670,2150,1650,5900
    im3xmin,im3xmax,im3ymin,im3ymax = 0,4920,5900,11999
    im4xmin,im4xmax,im4ymin,im4ymax = 2600,4920,0,5900
    im5xmin,im5xmax,im5ymin,im5ymax = 5870,10300,0,11999
    im6xmin,im6xmax,im6ymin,im6ymax = 11260,11930,0,11999

if glstart == 120:
    im1xmin,im1xmax,im1ymin,im1ymax = 330,1030,0,11999
    im2xmin,im2xmax,im2ymin,im2ymax = 1440,4630,8000,11999
    im3xmin,im3xmax,im3ymin,im3ymax = 1970,4630,4430,8000
    im4xmin,im4xmax,im4ymin,im4ymax = 1970,4250,0,4430
    im5xmin,im5xmax,im5ymin,im5ymax = 4630,11830,0,11999


##################################################
# Add random noise
##################################################

im1 = img[im1ymin:im1ymax,im1xmin:im1xmax]
im1r = im1 + (np.random.poisson(lam=2.,size=np.shape(im1))/200.)
fits.HDUList([fits.PrimaryHDU(im1r)]).writeto('../combmaps12000/'+str(glstart)+'_im1_rand_deadtime.fits')

if glstart == 20 or glstart == 60 or glstart == 80 or glstart == 100 or glstart == 120:
    im2 = img[im2ymin:im2ymax,im2xmin:im2xmax]
    im2r = im2 + (np.random.poisson(lam=2.,size=np.shape(im2))/200.)
    fits.HDUList([fits.PrimaryHDU(im2r)]).writeto('../combmaps12000/'+str(glstart)+'_im2_rand_deadtime.fits')
    im3 = img[im3ymin:im3ymax,im3xmin:im3xmax]
    im3r = im3 + (np.random.poisson(lam=2.,size=np.shape(im3))/200.)
    fits.HDUList([fits.PrimaryHDU(im3r)]).writeto('../combmaps12000/'+str(glstart)+'_im3_rand_deadtime.fits')
    im4 = img[im4ymin:im4ymax,im4xmin:im4xmax]
    im4r = im4 + (np.random.poisson(lam=2.,size=np.shape(im4))/200.)
    fits.HDUList([fits.PrimaryHDU(im4r)]).writeto('../combmaps12000/'+str(glstart)+'_im4_rand_deadtime.fits')
    im5 = img[im5ymin:im5ymax,im5xmin:im5xmax]
    im5r = im5 + (np.random.poisson(lam=2.,size=np.shape(im5))/200.)
    fits.HDUList([fits.PrimaryHDU(im5r)]).writeto('../combmaps12000/'+str(glstart)+'_im5_rand_deadtime.fits')

if glstart == 80:
    im6 = img[im6ymin:im6ymax,im6xmin:im6xmax]
    im6r = im6 + (np.random.poisson(lam=2.,size=np.shape(im6))/200.)
    fits.HDUList([fits.PrimaryHDU(im6r)]).writeto('../combmaps12000/'+str(glstart)+'_im6_rand_deadtime.fits')

print 'Added noise to all images'

##################################################
# Run sextractor
##################################################
os.system('sex ../combmaps12000/'+str(glstart)+'_im1_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im1_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im1_rand_deadtime.fits')

if glstart == 20 or glstart == 60 or glstart == 80 or glstart == 100 or glstart == 120:
    os.system('sex ../combmaps12000/'+str(glstart)+'_im2_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im2_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im2_rand_deadtime.fits')

    os.system('sex ../combmaps12000/'+str(glstart)+'_im3_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im3_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im3_rand_deadtime.fits')

    os.system('sex ../combmaps12000/'+str(glstart)+'_im4_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im4_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im4_rand_deadtime.fits')

    os.system('sex ../combmaps12000/'+str(glstart)+'_im5_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im5_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im5_rand_deadtime.fits')

if glstart == 80:
    os.system('sex ../combmaps12000/'+str(glstart)+'_im6_rand_deadtime.fits -c ~/sextractor/daofind.sex -CATALOG_NAME ../combmaps12000/sex_'+str(glstart)+'_im6_rand_deadtime.fits -CHECKIMAGE_NAME ../combmaps12000/background_'+str(glstart)+'_im6_rand_deadtime.fits')

##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################
'''
for a in range(0, 360, 20):
    b = a + 20
    try:
        data = fits.open('../combmaps12000/sex_'+str(a)+'to'+str(b)+'.fits')[1].data
    except IOError:
        print 'IOError'
        continue
    dgl = data.X_IMAGE * 20./12000. + a
    dgb = data.Y_IMAGE * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(a)+'to'+str(b)+'_edit.fits')
'''
# Or do this for one run

im1sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im1_rand_deadtime.fits')[1].data
data = im1sex
xfac = im1xmin
yfac = im1ymin
dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im1_rand_deadtime_edit.fits')

if glstart == 20 or glstart == 60 or glstart == 80 or glstart == 100 or glstart == 120:
    im2sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im2_rand_deadtime.fits')[1].data
    data = im2sex
    xfac = im2xmin
    yfac = im2ymin
    dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
    dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im2_rand_deadtime_edit.fits')

    im3sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im3_rand_deadtime.fits')[1].data
    data = im3sex
    xfac = im3xmin
    yfac = im3ymin
    dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
    dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im3_rand_deadtime_edit.fits')
    
    im4sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im4_rand_deadtime.fits')[1].data
    data = im4sex
    xfac = im4xmin
    yfac = im4ymin
    dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
    dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im4_rand_deadtime_edit.fits')

    im5sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im5_rand_deadtime.fits')[1].data
    data = im5sex
    xfac = im5xmin
    yfac = im5ymin
    dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
    dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im5_rand_deadtime_edit.fits')

if glstart == 80:
    im6sex = fits.open('../combmaps12000/sex_'+str(glstart)+'_im6_rand_deadtime.fits')[1].data
    data = im6sex
    xfac = 11260
    yfac = 0
    dgl = (data.X_IMAGE+xfac) * 20./12000. + glstart
    dgb = (data.Y_IMAGE+yfac) * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb',  format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('../combmaps12000/sex_'+str(glstart)+'_im6_rand_deadtime_edit.fits')

print 'Converted coordinates to gl/gb' 

##################################################
# Combine all tables
##################################################
#sexlist = np.loadtxt('../combmaps12000/sextractorlist.txt', dtype='str')
#tottable = Table.read('../combmaps12000/sex_0to20_edit.fits', format='fits')

tottable = Table.read('../combmaps12000/sex_'+str(glstart)+'_im1_rand_deadtime_edit.fits',format='fits')

if glstart == 40:
    sexlist = []

if glstart == 80:
    sexlist = ['sex_'+str(glstart)+'_im2_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im3_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im4_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im5_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im6_rand_deadtime_edit.fits']
if glstart == 20 or glstart == 60 or glstart == 100 or glstart == 120:
    sexlist = ['sex_'+str(glstart)+'_im2_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im3_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im4_rand_deadtime_edit.fits','sex_'+str(glstart)+'_im5_rand_deadtime_edit.fits']

for i in sexlist:
    if len(sexlist) > 0:
        tablen = Table.read('../combmaps12000/'+i, format='fits')
        tottable = vstack([tottable, tablen])

ascii.write(tottable, '../sex_total_'+str(glstart)+'_rand_deadtime.txt')

print 'Combined all data tables'

##################################################
# Now match with Tycho 2/2MASS catalogs
##################################################
sex = Table.read('../sex_total_'+str(glstart)+'_rand_deadtime.txt', format='ascii')
t2 = Table.read('../tycho2_2mass_matches.txt', format='ipac')

sexgal = SkyCoord(sex['gl']*u.degree, sex['gb']*u.degree, frame='galactic')
t2gal = SkyCoord(t2['gl_01']*u.degree, t2['gb_01']*u.degree, frame='galactic')

t2ind, sexind, angsep, dist3d = search_around_sky(t2gal, sexgal, 6.*u.arcsec)

plt.hist(angsep*3600., bins=100), plt.show()

sex = sex[sexind]
t3 = t2[t2ind]

print len(sex)

combtable = hstack([sex, t3])
combtable.rename_column('galex_id_01', 'galex_id')
combtable.rename_column('cat_id_01', 'cat_id')
combtable.rename_column('fuv_cps_01', 'fuv_cps')
combtable.rename_column('nuv_cps_01', 'nuv_cps')
combtable.rename_column('vjmag_01', 'VJmag')
combtable.rename_column('bjmag_01', 'BJmag')
combtable.rename_column('j_m', 'j')
combtable.rename_column('h_m', 'h')
combtable.rename_column('k_m', 'k')


ascii.write(combtable, '../sextractor_'+str(glstart)+'_rand_deadtime_t2_2mass.txt')

print 'Matched with T2 and 2MASS, finished'
print 'Total objects matched =', len(combtable)
