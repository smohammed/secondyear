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


#sg = fits.open('../sex-gaia.fits')[1].data
sg = Table.read('../sex_vphas_gaia.txt', format='ascii')
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

ebv = []

count = 0

for i in range(0, len(sg), 5000):
    print i
    sgcut = sg[i:i+5000]
    sggl = sgcut['gl_sex'].tolist()
    sggb = sgcut['gb_sex'].tolist()
    dmmatchindcut = dmmatchind[i:i+5000]

    gsf = query(sggl, sggb, coordsys='gal', mode='lite')

    for line in range(0, 5000):
        try:
            ebv.append(gsf['best'][line][dmmatchindcut[line]])
        except IndexError:
            print 'index error'
            count += 1
            pass

ebv = np.array(ebv)

sg = Table(sg)
sg['dist'] = pc
sg['distmod'] = distmod
sg['Mg'] = Mg
sg['DMdust'] = dmmatch
sg['ebv'] = ebv

print 'takes ~10 minutes to save...'
#ascii.write(sg, 'sex-gaia-dust.txt', format='ipac')
ascii.write(sg, 'sex_vphas_gaia_dust.txt', format='ipac')
