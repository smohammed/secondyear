from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

#gaialist = np.loadtxt('../../GAIA/planefiles/alist_planefiles.txt', dtype='str')
gaialist = np.loadtxt('../../GAIA/DR1/filelist.txt', dtype='str')

catalog = fits.open('../../galexscans/starcat_allscans_03-26-18.fits')[1].data
cut = np.where(catalog['expsum'] > 10.)
catalog = catalog[cut]

catgal = SkyCoord(catalog['gl']*u.deg, catalog['gb']*u.deg, frame='galactic')

table = Table()

for line in range(len(gaialist)):
	print gaialist[line]
	#gaia = Table.read('../../GAIA/planefiles/'+gaialist[line], format='ascii')
	gaia = fits.open('../../GAIA/DR1/'+gaialist[line])[1].data
	gaiagal = SkyCoord(gaia['l']*u.deg, gaia['b']*u.deg, frame='galactic')
	catind, gaiaind, angsep, ang3d = search_around_sky(catgal, gaiagal, 3*u.arcsec)
	if len(gaiaind) > 0:
		g2 = Table(gaia[gaiaind])
		g2['angsep'] = angsep
		table = vstack([table, g2])
		#ascii.write(g2, '../../GAIA/galexmatches/'+gaialist[line][11:22]+'_gaia_galexmatch.txt', format='basic')
		ascii.write(g2, '../../GAIA/galexmatches/'+gaialist[line][11:22]+'_gaia_galexmatch.txt', format='basic')
		print len(g2)

ascii.write(table, 'gaia_sextractor_match.txt', format='basic')