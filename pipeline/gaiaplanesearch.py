from astropy.io import ascii
import numpy as np
from astropy.table import Table, vstack

gaialist = np.loadtxt('../../GAIA/DR1/filelist_for_planesearch.txt', dtype='str')

files = []

for line in range(len(gaialist)):
	print gaialist[line]
	gaia = Table.read('../../GAIA/DR1/'+gaialist[line])
	planecut = np.where(np.abs(gaia['b']) < 10.)
	if len(planecut[0] > 0.):
		files.append(gaialist[line])
		gaia = gaia[planecut]
		ascii.write(gaia, '../../GAIA/planefiles/'+gaialist[line]+'planefiles.txt', format='basic')

ascii.write(files, 'gaiadr1planecutlist.txt', format='basic')