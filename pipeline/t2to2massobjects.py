from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs7


tb = Table.read('../tychobstarmatch.fits', format='fits')
tbcut = np.where((tb['gb'] > -10) & (tb['gb'] < 10))
tb2 = tb[tbcut]

nuvlowlim = 1000
nuvhilim = 10000

tbnuvcut2 = np.where((tb2['nuv_cps'] > nuvlowlim) & (tb2['nuv_cps'] < nuvhilim))

tb3 = tb2[tbnuvcut2]

print 'now writing file'
ascii.write(tb3, '../tycho2_nuvcut_'+str(nuvlowlim)+'-'+str(nuvhilim)+'_gbcut.txt', format='ipac')

ra = tb3['ra']
dec = tb3['dec']

newtable = Table([ra, dec], names=('ra', 'dec'))

ascii.write(newtable, '../tycho2_nuvcut_'+str(nuvlowlim)+'-'+str(nuvhilim)+'_gbcut_posonly.txt', format='ipac')
