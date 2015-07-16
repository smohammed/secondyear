from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20


##################################################
# Get output from sextractor, convert to gl,gb,nuv
##################################################

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


##################################################
# Combine all tables
##################################################

sexlist = np.loadtxt('../combmaps12000/sextractorlist.txt', dtype='str')
tottable = Table.read('../combmaps12000/sex_0to20_edit.fits', format='fits')

for i in sexlist:
    tablen = Table.read('../combmaps12000/'+i, format='fits')
    tottable = vstack([tottable, tablen])

ascii.write(tottable, '../sex_total.txt')

##################################################
# Now match with Tycho 2/2MASS catalogs
##################################################
sex = Table.read('../sex_total.txt', format='ascii')
t2 = Table.read('../tycho2_2mass_matches.txt', format='ipac')

sexgal = SkyCoord(sex['gl']*u.deg, sex['gb']*u.deg, frame='galactic')
t2gal = SkyCoord(t2['gl_01']*u.degree, t2['gb_01']*u.degree, frame='galactic')

t2ind, sexind, angsep, dist3d = search_around_sky(t2gal, sexgal, 6.*u.arcsec)

plt.hist(angsep*3600., bins=100), plt.show()

sex = sex[sexind]
t2 = t2[t2ind]

combtable = hstack([sex, t2])
ascii.write(combtable, '../sextractor_t2_2mass.txt')

plt.imshow(intmap,vmin=0,vmax=0.7,origin='lower',extent=[120,140,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray)