from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

rarange = '275-300'

# Assign variables
ps1 = fits.open('../ps1/ps1_ra'+rarange+'_g10-20_unedited_redo.fit')[1].data
ps1 = Table(ps1)
gra = ps1['gra']
rra = ps1['rra']
ira = ps1['ira']
zra = ps1['zra']
yra = ps1['yra']

gdec = ps1['gdec']
rdec = ps1['rdec']
idec = ps1['idec']
zdec = ps1['zdec']
ydec = ps1['ydec']

g_ps = ps1['gPSFMag']
r_ps = ps1['rPSFMag']
i_ps = ps1['iPSFMag']
z_ps = ps1['zPSFMag']
y_ps = ps1['yPSFMag']

# Remove all null data
x = np.where((gra == -999.) | (rra == -999.) | (ira == -999.) | (zra == -999.) | (yra == -999.) | (gdec == -999.) | (rdec == -999.) | (idec == -999.) | (zdec == -999.) | (ydec == -999.) | (g_ps == -999.) | (r_ps == -999.) | (i_ps == -999.) | (z_ps == -999.) | (y_ps == -999.))
ps1.remove_rows(x)

# Calculate mean of ra/dec from band coordinates
ramean = (ps1['gra'] + ps1['rra'] + ps1['ira'] + ps1['zra'] + ps1['yra'])/5.
decmean = (ps1['gdec'] + ps1['rdec'] + ps1['idec'] + ps1['zdec'] + ps1['ydec'])/5.

psgal = SkyCoord(ramean*u.deg, decmean*u.deg, frame='icrs').galactic

# Make table
cols = fits.ColDefs([fits.Column(name='ramean', format='D', array=ramean), fits.Column(name='decmean', format='D', array=decmean), fits.Column(name='glmean', format='D', array=psgal.l.degree), fits.Column(name='gbmean', format='D', array=psgal.b.degree), fits.Column(name='g_ps', format='D', array=ps1['gPSFMag']), fits.Column(name='r_ps', format='D', array=ps1['rPSFMag']), fits.Column(name='i_ps', format='D', array=ps1['iPSFMag']), fits.Column(name='z_ps', format='D', array=ps1['zPSFMag']), fits.Column(name='y_ps', format='D', array=ps1['yPSFMag'])])

print 'making table'

endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('../ps1/ps1_ra'+rarange+'_g10-20.fits')
