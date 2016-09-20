from astropy.io import fits, ascii
from astropy.table import Table, hstack
import numpy as np
from dustquery import query
import math

def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


sg = fits.open('../sex-gaia.fits')[1].data
pc = 1000 / sg['parallax']
negpar = np.where(pc > 0)
pc = pc[negpar]
sg = sg[negpar]
distmod = 5 * np.log10(pc) - 5
Mg = sg['phot_g_mean_mag'] - distmod

DM_dust = np.arange(4, 19.1, 0.5)
dmmatch, dmmatchind = [], []

for i in range(len(distmod)):
    match = find_nearest(DM_dust, distmod[i])
    dmmatch.append(match)
    dmmatchind.append(np.where(DM_dust == match)[0][0])

dmmatch = np.array(dmmatch)
dmmatchind = np.array(dmmatchind)

sggl = sg['gl'].tolist()
sggb = sg['gb'].tolist()

ebv = []

for i in range(len(sg))[:1000]:
    gsf = query(sggl[i], sggb[i], coordsys='gal', mode='lite')
    ebv.append(gsf['best'][dmmatchind[i]])

ebv = np.array(ebv)
ebvfinal = Table([ebv])

ebvfinal.rename_column('col0', 'ebv')
ascii.write(ebvfinal, 'ebv_sex_gaia.txt', format='basic')



#[u'b', u'GR',u'success',u'l',u'best',u'DM_reliable_max',u'ra',u'samples',u'n_stars',u'converged',u'dec',u'DM_reliable_min',u'distmod']


