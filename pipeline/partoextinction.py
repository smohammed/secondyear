from astropy.io import fits
from astropy.table import Table
import numpy as np
from dustquery import query
import math


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


#sg = fits.open('../sex_gaia_fec.fits')[1].data
sg = fits.open('../sex_gaia_allnuv.fits')[1].data
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

#sg.rename_column('ra', 'ra_gaia')
#sg.rename_column('dec', 'dec_gaia')

print 'takes ~10 minutes to save...'
#ascii.write(sg, 'sex-gaia-dust.txt', format='ipac')
#ascii.write(sg, 'sex_vphas_gaia_dust.txt', format='ipac')

#fits.Column(name='fec', format='20A', array=sg['fec'])

endtable = fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='NUMBER', format='E', array=sg['NUMBER']), fits.Column(name='X_IMAGE', format='E', array=sg['X_IMAGE']), fits.Column(name='Y_IMAGE', format='E', array=sg['Y_IMAGE']), fits.Column(name='FLUX_APER', format='E', array=sg['FLUX_APER']), fits.Column(name='A_IMAGE', format='E', array=sg['A_IMAGE']), fits.Column(name='B_IMAGE', format='E', array=sg['B_IMAGE']), fits.Column(name='THETA_IMAGE', format='E', array=sg['THETA_IMAGE']), fits.Column(name='FWHM_IMAGE', format='E', array=sg['FWHM_IMAGE']), fits.Column(name='nuv_mag', format='E', array=sg['nuv_mag']), fits.Column(name='gl_sex', format='E', array=sg['gl_sex']), fits.Column(name='gb_sex', format='E', array=sg['gb_sex']), fits.Column(name='ra_sex', format='E', array=sg['ra_1']), fits.Column(name='dec_sex', format='E', array=sg['dec_sex']), fits.Column(name='hip', format='J', array=sg['hip']), fits.Column(name='tycho2_id', format='12A', array=sg['tycho2_id']), fits.Column(name='solution_id', format='K', array=sg['solution_id']), fits.Column(name='source_id', format='K', array=sg['source_id']), fits.Column(name='random_index', format='K', array=sg['random_index']), fits.Column(name='ref_epoch', format='D', array=sg['ref_epoch']), fits.Column(name='ra_gaia', format='D', array=sg['ra_gaia']), fits.Column(name='ra_error', format='D', array=sg['ra_error']), fits.Column(name='dec_gaia', format='D', array=sg['dec_gaia']), fits.Column(name='dec_error', format='D', array=sg['dec_error']), fits.Column(name='parallax', format='D', array=sg['parallax']), fits.Column(name='parallax_error', format='D', array=sg['parallax_error']), fits.Column(name='pmra', format='D', array=sg['pmra']), fits.Column(name='pmra_error', format='D', array=sg['pmra_error']), fits.Column(name='pmdec', format='D', array=sg['pmdec']), fits.Column(name='pmdec_error', format='D', array=sg['pmdec_error']), fits.Column(name='ra_dec_corr', format='E', array=sg['ra_dec_corr']), fits.Column(name='ra_parallax_corr', format='E', array=sg['ra_parallax_corr']), fits.Column(name='ra_pmra_corr', format='E', array=sg['ra_pmra_corr']), fits.Column(name='ra_pmdec_corr', format='E', array=sg['ra_pmdec_corr']), fits.Column(name='dec_parallax_corr', format='E', array=sg['dec_parallax_corr']), fits.Column(name='dec_pmra_corr', format='E', array=sg['dec_pmra_corr']), fits.Column(name='dec_pmdec_corr', format='E', array=sg['dec_pmdec_corr']), fits.Column(name='parallax_pmra_corr', format='E', array=sg['parallax_pmra_corr']), fits.Column(name='parallax_pmdec_corr', format='E', array=sg['parallax_pmdec_corr']), fits.Column(name='pmra_pmdec_corr', format='E', array=sg['pmra_pmdec_corr']), fits.Column(name='astrometric_n_obs_al', format='J', array=sg['astrometric_n_obs_al']), fits.Column(name='astrometric_n_obs_ac', format='J', array=sg['astrometric_n_obs_ac']), fits.Column(name='astrometric_n_good_obs_al', format='J', array=sg['astrometric_n_good_obs_al']), fits.Column(name='astrometric_n_good_obs_ac', format='J', array=sg['astrometric_n_good_obs_ac']), fits.Column(name='astrometric_n_bad_obs_al', format='J', array=sg['astrometric_n_bad_obs_al']), fits.Column(name='astrometric_n_bad_obs_ac', format='J', array=sg['astrometric_n_bad_obs_ac']), fits.Column(name='astrometric_delta_q', format='E', array=sg['astrometric_delta_q']), fits.Column(name='astrometric_excess_noise', format='D', array=sg['astrometric_excess_noise']), fits.Column(name='astrometric_excess_noise_sig', format='D', array=sg['astrometric_excess_noise']), fits.Column(name='astrometric_primary_flag', format='L', array=sg['astrometric_primary_flag']), fits.Column(name='astrometric_relegation_factor', format='E', array=sg['astrometric_relegation_factor']), fits.Column(name='astrometric_weight_al', format='E', array=sg['astrometric_weight_al']), fits.Column(name='astrometric_weight_ac', format='E', array=sg['astrometric_weight_ac']), fits.Column(name='astrometric_priors_used', format='J', array=sg['astrometric_priors_used']), fits.Column(name='matched_observations', format='I', array=sg['matched_observations']), fits.Column(name='duplicated_source', format='L', array=sg['duplicated_source']), fits.Column(name='scan_direction_strength_k1', format='E', array=sg['scan_direction_strength_k1']), fits.Column(name='scan_direction_strength_k2', format='E', array=sg['scan_direction_strength_k2']), fits.Column(name='scan_direction_strength_k3', format='E', array=sg['scan_direction_strength_k3']), fits.Column(name='scan_direction_strength_k4', format='E', array=sg['scan_direction_strength_k4']), fits.Column(name='scan_direction_mean_k1', format='E', array=sg['scan_direction_mean_k1']), fits.Column(name='scan_direction_mean_k2', format='E', array=sg['scan_direction_mean_k2']), fits.Column(name='scan_direction_mean_k3', format='E', array=sg['scan_direction_mean_k3']), fits.Column(name='scan_direction_mean_k4', format='E', array=sg['scan_direction_mean_k4']), fits.Column(name='phot_g_n_obs', format='J', array=sg['phot_g_n_obs']), fits.Column(name='phot_g_mean_flux', format='D', array=sg['phot_g_mean_flux']), fits.Column(name='phot_g_mean_flux_error', format='D', array=sg['phot_g_mean_flux_error']), fits.Column(name='phot_g_mean_mag', format='D', array=sg['phot_g_mean_mag']), fits.Column(name='phot_variable_flag', format='13A', array=sg['phot_variable_flag']), fits.Column(name='gl_gaia', format='D', array=sg['gl_gaia']), fits.Column(name='gb_gaia', format='D', array=sg['gb_gaia']), fits.Column(name='ecl_lon', format='D', array=sg['ecl_lon']), fits.Column(name='ecl_lat', format='D', array=sg['ecl_lat']), fits.Column(name='angsep_sg', format='E', array=sg['angsep_sg']), fits.Column(name='dist', format='E', array=pc), fits.Column(name='distmod', format='E', array=distmod), fits.Column(name='Mg', format='E', array=Mg), fits.Column(name='DMdust', format='E', array=dmmatch), fits.Column(name='ebv', format='E', array=ebv)]))
endtable.writeto('../sex_gaia_dust_allnuv.fits')
