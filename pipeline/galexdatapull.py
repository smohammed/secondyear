import numpy as np
from astropy.io import ascii
from astropy.table import Table, hstack, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

galexlist = np.loadtxt('../filelists/galexfilelist.txt', dtype='string')
path = '/home/cal/smohammed/Downloads/'
data = Table.read(path+galexlist[0], format='fits')

ggoid_dec = Table([data['ggoid_dec']])
alpha_j2000 = Table([data['alpha_j2000']])
delta_j2000 = Table([data['delta_j2000']])
E_bv = Table([data['E_bv']])
glon = Table([data['glon']])
glat = Table([data['glat']])
nuv = Table([data['nuv_mag']])
fuv = Table([data['fuv_mag']])
NUV_FLUX_AUTO = Table([data['NUV_FLUX_AUTO']])
NUV_FLUXERR_AUTO = Table([data['NUV_FLUXERR_AUTO']])
FUV_FLUX_AUTO = Table([data['FUV_FLUX_AUTO']])
FUV_FLUXERR_AUTO = Table([data['FUV_FLUXERR_AUTO']])

for i in range(1, len(galexlist), 1):
    tablen = Table.read(path+galexlist[i], format='fits')
    ggoid_dec = vstack([ggoid_dec, Table([tablen['ggoid_dec']])])
    alpha_j2000 = vstack([alpha_j2000, Table([tablen['alpha_j2000']])])
    delta_j2000 = vstack([delta_j2000, Table([tablen['delta_j2000']])])
    E_bv = vstack([E_bv, Table([tablen['E_bv']])])
    glon = vstack([glon, Table([tablen['glon']])])
    glat = vstack([glat, Table([tablen['glat']])])
    nuv = vstack([nuv, Table([tablen['nuv_mag']])])
    fuv = vstack([fuv, Table([tablen['fuv_mag']])])
    NUV_FLUX_AUTO = vstack([NUV_FLUX_AUTO, Table([tablen['NUV_FLUX_AUTO']])])
    NUV_FLUXERR_AUTO = vstack([NUV_FLUXERR_AUTO, Table([tablen['NUV_FLUXERR_AUTO']])])
    FUV_FLUX_AUTO = vstack([FUV_FLUX_AUTO, Table([tablen['FUV_FLUX_AUTO']])])
    FUV_FLUXERR_AUTO = vstack([FUV_FLUXERR_AUTO, Table([tablen['FUV_FLUXERR_AUTO']])])
    print i

finaltable = hstack([ggoid_dec, alpha_j2000, delta_j2000, E_bv, glon, glat, nuv, fuv, NUV_FLUX_AUTO, NUV_FLUXERR_AUTO, FUV_FLUX_AUTO, FUV_FLUXERR_AUTO])

ascii.write(finaltable, '../galex120data.txt')

t2 = Table.read('../tycho2_2mass_matches.txt', format='ascii')
galexgal = SkyCoord(finaltable['alpha_j2000']*u.deg, finaltable['delta_j2000']*u.deg, frame='icrs')
t2gal = SkyCoord(t2['ra_t2']*u.deg, t2['dec_t2']*u.deg, frame='icrs')

t2ind, galind, angsep, dist3d = search_around_sky(t2gal, galexgal, 1.*u.arcsec)

t3 = t2[t2ind]
galex2 = finaltable[galind]

alldata = hstack([galex2, t3])

alldata.rename_column('glon', 'gl')
alldata.rename_column('glat', 'gb')

ascii.write(alldata, '../galex120_2mass_t2.txt')
