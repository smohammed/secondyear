from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

# So far only added ra 0-200
#tables = ['ps1_ra0-25_g10-20.fits', 'ps1_ra25-50_g10-20.fits', 'ps1_ra50-75_g10-20.fits', 'ps1_ra75-100_g10-20.fits']

#tables = ['ps1_ra100-125_g10-20.fits', 'ps1_ra125-150_g10-20.fits', 'ps1_ra150-175_g10-20.fits', 'ps1_ra175-200_g10-20.fits']

tables = ['ps1_ra200-225_g10-20.fits', 'ps1_ra225-250_g10-20.fits', 'ps1_ra250-275_g10-20.fits', 'ps1_ra275-300_g10-20.fits']

#tables = ['ps1_ra300-325_g10-20.fits', 'ps1_ra325-350_g10-20.fits', 'ps1_ra350-360_g10-20.fits']

glrange = '200-300'

comb = Table()

for tbl in tables:
    print tbl
    table = fits.open('../ps1/'+tbl)[1].data
    cut = np.where((table['gbmean'] > -10) & (table['gbmean'] < 10))
    table = Table(table[cut])
    comb = vstack([comb, table])

cols = fits.ColDefs([fits.Column(name='ramean', format='D', array=comb['ramean']), fits.Column(name='decmean', format='D', array=comb['decmean']), fits.Column(name='glmean', format='D', array=comb['glmean']), fits.Column(name='gbmean', format='D', array=comb['gbmean']), fits.Column(name='g_ps', format='D', array=comb['g_ps']), fits.Column(name='r_ps', format='D', array=comb['r_ps']), fits.Column(name='i_ps', format='D', array=comb['i_ps']), fits.Column(name='z_ps', format='D', array=comb['z_ps']), fits.Column(name='y_ps', format='D', array=comb['y_ps'])])

print 'making PS table'

endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('../ps1/ps1_planecut_'+glrange+'_g10-20.fits')


ps = fits.open('../ps1/ps1_planecut_'+glrange+'_g10-20.fits')[1].data
psgal = SkyCoord(ps['glmean']*u.deg, ps['gbmean']*u.deg, frame='galactic')

cat = Table(fits.open('../starcat_fwhm_12-14.fits')[1].data)
catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')

catind, psind, angsep, ang3d = search_around_sky(catgal, psgal, 3*u.arcsec)
c2 = cat[catind]
ps2 = ps[psind]
comb = hstack([Table(c2), Table(ps2)])

comb['angsep'] = angsep
comb.rename_column('gl', 'gl_sex')
comb.rename_column('gb', 'gb_sex')
comb.rename_column('ra', 'ra_sex')
comb.rename_column('dec', 'dec_sex')
comb.rename_column('ramean', 'ra_ps')
comb.rename_column('decmean', 'dec_ps')
comb.rename_column('glmean', 'gl_ps')
comb.rename_column('gbmean', 'gb_ps')


cols = fits.ColDefs([fits.Column(name='NUMBER', format='D', array=comb['NUMBER']), fits.Column(name='X_IMAGE', format='D', array=comb['X_IMAGE']), fits.Column(name='Y_IMAGE', format='D', array=comb['Y_IMAGE']), fits.Column(name='FLUX_AUTO', format='D', array=comb['FLUX_AUTO']), fits.Column(name='FLUXERR_AUTO', format='D', array=comb['FLUXERR_AUTO']), fits.Column(name='FLUX_APER', format='D', array=comb['FLUX_APER']), fits.Column(name='A_IMAGE', format='D', array=comb['A_IMAGE']), fits.Column(name='B_IMAGE', format='D', array=comb['B_IMAGE']), fits.Column(name='THETA_IMAGE', format='D', array=comb['THETA_IMAGE']), fits.Column(name='FWHM_IMAGE', format='D', array=comb['FWHM_IMAGE']), fits.Column(name='x_new', format='D', array=comb['x_new']), fits.Column(name='y_new', format='D', array=comb['y_new']), fits.Column(name='nuv', format='D', array=comb['nuv']), fits.Column(name='gl_sex', format='D', array=comb['gl_sex']), fits.Column(name='gb_sex', format='D', array=comb['gb_sex']), fits.Column(name='ra_sex', format='D', array=comb['ra_sex']), fits.Column(name='dec_sex', format='D', array=comb['dec_sex']), fits.Column(name='ra_ps', format='D', array=comb['ra_ps']), fits.Column(name='dec_ps', format='D', array=comb['dec_ps']), fits.Column(name='gl_ps', format='D', array=comb['gl_ps']), fits.Column(name='gb_ps', format='D', array=comb['gb_ps']), fits.Column(name='g_ps', format='D', array=comb['g_ps']), fits.Column(name='r_ps', format='D', array=comb['r_ps']), fits.Column(name='i_ps', format='D', array=comb['i_ps']), fits.Column(name='z_ps', format='D', array=comb['z_ps']), fits.Column(name='y_ps', format='D', array=comb['y_ps']), fits.Column(name='angsep_sps', format='D', array=comb['angsep'])])

print 'making sextractor combined table'

endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('../ps1/sex_ps1_'+glrange+'_g10-20.fits')
