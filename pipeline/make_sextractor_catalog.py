from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy.table import vstack
from astropy.io import ascii

for a in range(0, 360, 20):
    b = a + 20
    try:
        data = fits.open('sex_'+str(a)+'to'+str(b)+'.fits')[1].data
    except IOError:
        continue
    dgl = data.X_IMAGE * 20./12000. + a
    dgb = data.Y_IMAGE * 20./12000. - 10.
    nuv = -2.5*np.log10(data.FLUX_AUTO*10.) + 20.08
    new_cols = fits.ColDefs([fits.Column(name='gl', format='1E', array=dgl), fits.Column(name='gb', format='1E', array=dgb), fits.Column(name='nuv', format='1E', array=nuv)])
    hdu = fits.BinTableHDU.from_columns(data.columns + new_cols)
    hdu.writeto('sex_'+str(a)+'to'+str(b)+'_edit.fits')

sexlist = np.loadtxt('sextractorlist.txt', dtype='str')
tottable = Table.read('sex_0to20_edit.fits', format='fits')

for i in sexlist:
    tablen = Table.read(i, format='fits')
    tottable = vstack([tottable, tablen])

ascii.write(tottable, 'sex_total.txt')

sex = Table.read('sex_total.txt', format='ascii')
t2 = Table.read('tycho2_2mass_matches.txt', format='ipac')

sexgal = SkyCoord(sex['gl']*u.deg, sex['gb']*u.deg, frame='galactic')
t2gal = SkyCoord(t2['gl_01']*u.degree, t2['gb_01']*u.degree, frame='galactic')

t2ind, sexind, angsep, dist3d = search_around_sky(sexgal, t2gal, 6.*u.arcsec)

plt.hist(angsep*3600., bins=100), plt.show()

