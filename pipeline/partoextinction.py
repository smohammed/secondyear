from astropy.io import fits
from astropy.table import Table
import numpy as np
from dustquery import query
import math
import scipy.interpolate as inp


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

cat = 1
gais = 0

#################################################################
# For SExtractor catalog
#################################################################
if cat == 1:
    sg = Table.read('../sex_tgas_comb_02-28-2018.txt', format='ascii')
    #sg = fits.open('../galexplane_tgas.fits')[1].data

    # Calculate parameters to get distance
    oldpar = sg['parallax']
    parallax = oldpar * (0.5 + 0.5 * np.sqrt(1. - 16 * (sg['parallax_error']**2/oldpar**2)))
    pc = 1000. / parallax
    negpar = np.where(pc > 0)
    pc = pc[negpar]
    parallax = parallax[negpar]
    sg = sg[negpar]
    distmod = 5. * np.log10(pc) - 5.
    Mg = sg['phot_g_mean_mag'] - distmod
    DM_dust = np.arange(4, 19.1, 0.5)

    # Use every 5000 iterations because that's the max you can query
    ebv = []
    count = 0
    for i in range(0, len(sg), 5000):
        print(i)
        sgcut = sg[i:i+5000]
        sggl = sgcut['gl'].tolist()
        sggb = sgcut['gb'].tolist()
        distmodcut = distmod[i:i+5000]

        gsf = query(sggl, sggb, coordsys='gal', mode='lite')

        # Now match the distance indices to each ebv
        for line in range(0, 5000):
            try:
                dustfn = inp.interp1d(DM_dust, gsf['best'][line])
                ebv.append(dustfn(distmodcut[line]))

            except ValueError:
                count += 1
                ebv.append(0)

            except IndexError:
                print('index error')
                pass

    ebv = np.array(ebv)

    #cols = fits.ColDefs([fits.Column(name='X_IMAGE', format='D', array=sg['X_IMAGE']), fits.Column(name='Y_IMAGE', format='D', array=sg['Y_IMAGE']), fits.Column(name='FLUX_APER', format='D', array=sg['FLUX_APER']), fits.Column(name='A_IMAGE', format='D', array=sg['A_IMAGE']), fits.Column(name='B_IMAGE', format='D', array=sg['B_IMAGE']), fits.Column(name='THETA_IMAGE', format='D', array=sg['THETA_IMAGE']), fits.Column(name='FWHM_IMAGE', format='D', array=sg['FWHM_IMAGE']), fits.Column(name='nuv_mag', format='D', array=sg['nuv_mag']), fits.Column(name='gl_sex', format='D', array=sg['gl_sex']), fits.Column(name='gb_sex', format='D', array=sg['gb_sex']), fits.Column(name='ra_sex', format='D', array=sg['ra_sex']), fits.Column(name='dec_sex', format='D', array=sg['dec_sex']), fits.Column(name='hip', format='K', array=sg['hip']), fits.Column(name='tycho2_id', format='13A', array=sg['tycho2_id']), fits.Column(name='solution_id', format='K', array=sg['solution_id']), fits.Column(name='source_id', format='K', array=sg['source_id']), fits.Column(name='random_index', format='K', array=sg['random_index']), fits.Column(name='ref_epoch', format='D', array=sg['ref_epoch']), fits.Column(name='ra_gaia', format='D', array=sg['ra_gaia']), fits.Column(name='ra_error', format='D', array=sg['ra_error']), fits.Column(name='dec_gaia', format='D', array=sg['dec_gaia']), fits.Column(name='dec_error', format='D', array=sg['dec_error']), fits.Column(name='parallax', format='D', array=sg['parallax']), fits.Column(name='parallax_error', format='D', array=sg['parallax_error']), fits.Column(name='pmra', format='D', array=sg['pmra']), fits.Column(name='pmra_error', format='D', array=sg['pmra_error']), fits.Column(name='pmdec', format='D', array=sg['pmdec']), fits.Column(name='pmdec_error', format='D', array=sg['pmdec_error']), fits.Column(name='ra_dec_corr', format='E', array=sg['ra_dec_corr']), fits.Column(name='ra_parallax_corr', format='D', array=sg['ra_parallax_corr']), fits.Column(name='ra_pmra_corr', format='D', array=sg['ra_pmra_corr']), fits.Column(name='ra_pmdec_corr', format='D', array=sg['ra_pmdec_corr']), fits.Column(name='dec_parallax_corr', format='D', array=sg['dec_parallax_corr']), fits.Column(name='dec_pmra_corr', format='D', array=sg['dec_pmra_corr']), fits.Column(name='dec_pmdec_corr', format='D', array=sg['dec_pmdec_corr']), fits.Column(name='parallax_pmra_corr', format='D', array=sg['parallax_pmra_corr']), fits.Column(name='parallax_pmdec_corr', format='D', array=sg['parallax_pmdec_corr']), fits.Column(name='pmra_pmdec_corr', format='D', array=sg['pmra_pmdec_corr']), fits.Column(name='astrometric_n_obs_al', format='K', array=sg['astrometric_n_obs_al']), fits.Column(name='astrometric_n_obs_ac', format='K', array=sg['astrometric_n_obs_ac']), fits.Column(name='astrometric_n_good_obs_al', format='K', array=sg['astrometric_n_good_obs_al']), fits.Column(name='astrometric_n_good_obs_ac', format='K', array=sg['astrometric_n_good_obs_ac']), fits.Column(name='astrometric_n_bad_obs_al', format='K', array=sg['astrometric_n_bad_obs_al']), fits.Column(name='astrometric_n_bad_obs_ac', format='K', array=sg['astrometric_n_bad_obs_ac']), fits.Column(name='astrometric_delta_q', format='D', array=sg['astrometric_delta_q']), fits.Column(name='astrometric_excess_noise', format='D', array=sg['astrometric_excess_noise']), fits.Column(name='astrometric_excess_noise_sig', format='D', array=sg['astrometric_excess_noise']), fits.Column(name='astrometric_primary_flag', format='25A', array=sg['astrometric_primary_flag']), fits.Column(name='astrometric_relegation_factor', format='D', array=sg['astrometric_relegation_factor']), fits.Column(name='astrometric_weight_al', format='D', array=sg['astrometric_weight_al']), fits.Column(name='astrometric_weight_ac', format='D', array=sg['astrometric_weight_ac']), fits.Column(name='astrometric_priors_used', format='K', array=sg['astrometric_priors_used']), fits.Column(name='matched_observations', format='K', array=sg['matched_observations']), fits.Column(name='duplicated_source', format='18A', array=sg['duplicated_source']), fits.Column(name='scan_direction_strength_k1', format='D', array=sg['scan_direction_strength_k1']), fits.Column(name='scan_direction_strength_k2', format='D', array=sg['scan_direction_strength_k2']), fits.Column(name='scan_direction_strength_k3', format='D', array=sg['scan_direction_strength_k3']), fits.Column(name='scan_direction_strength_k4', format='D', array=sg['scan_direction_strength_k4']), fits.Column(name='scan_direction_mean_k1', format='D', array=sg['scan_direction_mean_k1']), fits.Column(name='scan_direction_mean_k2', format='D', array=sg['scan_direction_mean_k2']), fits.Column(name='scan_direction_mean_k3', format='D', array=sg['scan_direction_mean_k3']), fits.Column(name='scan_direction_mean_k4', format='D', array=sg['scan_direction_mean_k4']), fits.Column(name='phot_g_n_obs', format='J', array=sg['phot_g_n_obs']), fits.Column(name='phot_g_mean_flux', format='D', array=sg['phot_g_mean_flux']), fits.Column(name='phot_g_mean_flux_error', format='D', array=sg['phot_g_mean_flux_error']), fits.Column(name='phot_g_mean_mag', format='D', array=sg['phot_g_mean_mag']), fits.Column(name='phot_variable_flag', format='19A', array=sg['phot_variable_flag']), fits.Column(name='gl_gaia', format='D', array=sg['gl_gaia']), fits.Column(name='gb_gaia', format='D', array=sg['gb_gaia']), fits.Column(name='ecl_lon', format='D', array=sg['ecl_lon']), fits.Column(name='ecl_lat', format='D', array=sg['ecl_lat']), fits.Column(name='angsep_sg', format='D', array=sg['angsep_sg']), fits.Column(name='dist', format='D', array=pc), fits.Column(name='distmod', format='D', array=distmod), fits.Column(name='Mg', format='D', array=Mg), fits.Column(name='ebv', format='D', array=ebv), fits.Column(name='parallax_hogg', format='D', array=parallax)])

    cols = fits.ColDefs([fits.Column(name='nuv_mag', format='D', array=sg['nuv']), fits.Column(name='gl_sex', format='D', array=sg['gl']), fits.Column(name='gb_sex', format='D', array=sg['gb']), fits.Column(name='ra_sex', format='D', array=sg['ALPHA_J2000']), fits.Column(name='dec_sex', format='D', array=sg['DELTA_J2000']), fits.Column(name='phot_g_mean_mag', format='D', array=sg['phot_g_mean_mag']), fits.Column(name='gl_gaia', format='D', array=sg['l']), fits.Column(name='gb_gaia', format='D', array=sg['b']), fits.Column(name='parallax', format='D', array=sg['parallax']), fits.Column(name='parallax_error', format='D', array=sg['parallax_error']), fits.Column(name='dist', format='D', array=pc), fits.Column(name='distmod', format='D', array=distmod), fits.Column(name='MG', format='D', array=Mg), fits.Column(name='ebv', format='D', array=ebv), fits.Column(name='parallax_hogg', format='D', array=parallax)])

    endtable = fits.BinTableHDU.from_columns(cols)
    endtable.writeto('../galexplane_gaia_dust_interp_02-28-2018.fits')

#################################################################
# For GAIS instead
#################################################################
elif gais == 1:
    #sg = fits.open('../gais_gaia.fits')[1].data
    sg = fits.open('../gais_tgas.fits')[1].data
    sg = Table(sg)

    '''
    #sg.rename_column('Gaia_Parallax', 'parallax')
    #sg.rename_column('Gaia_Parallax_Err', 'parallax_error')
    #sg.rename_column('Gaia_G_Mag', 'phot_g_mean_mag')
    sg.rename_column('ra_1', 'ra_gais')
    sg.rename_column('dec_1', 'dec_gais')
    sg.rename_column('glon', 'gl_gais')
    sg.rename_column('glat', 'gb_gais')
    sg.rename_column('ra_2', 'ra_tgas')
    sg.rename_column('dec_2', 'dec_tgas')
    sg.rename_column('l', 'gl_tgas')
    sg.rename_column('b', 'gb_tgas')
    '''

    oldpar = sg['parallax']
    parallax = oldpar * (0.5 + 0.5 * np.sqrt(1. - 16 * (sg['parallax_error']**2/oldpar**2)))
    pc = 1000. / parallax
    negpar = np.where(pc > 0)
    pc = pc[negpar]
    parallax = parallax[negpar]
    sg = sg[negpar]
    distmod = 5. * np.log10(pc) - 5.
    Mg = sg['phot_g_mean_mag'] - distmod
    DM_dust = np.arange(4, 19.1, 0.5)

    # Use every 5000 iterations because that's the max you can query
    ebv = []
    count = 0
    for i in range(0, len(sg), 5000):
        print(i)
        sgcut = sg[i:i+5000]
        sggl = sgcut['gl_gais'].tolist()
        sggb = sgcut['gb_gais'].tolist()
        distmodcut = distmod[i:i+5000]

        gsf = query(sggl, sggb, coordsys='gal', mode='lite')

        # Now match the distance indices to each ebv
        for line in range(0, 5000):
            try:
                dustfn = inp.interp1d(DM_dust, gsf['best'][line])
                ebv.append(dustfn(distmodcut[line]))

            except ValueError:
                count += 1
                ebv.append(0)

            except IndexError:
                print('index error')
                pass

    ebv = np.array(ebv)

    cols = fits.ColDefs([fits.Column(name='nuv_mag', format='D', array=sg['nuv_mag']), fits.Column(name='nuv_magerr', format='D', array=sg['nuv_magerr']), fits.Column(name='fuv_mag', format='D', array=sg['fuv_mag']), fits.Column(name='fuv_magerr', format='D', array=sg['fuv_magerr']), fits.Column(name='gl_gais', format='D', array=sg['gl_gais']), fits.Column(name='gb_gais', format='D', array=sg['gb_gais']), fits.Column(name='ra_gais', format='D', array=sg['ra_gais']), fits.Column(name='dec_gais', format='D', array=sg['dec_gais']), fits.Column(name='hip', format='K', array=sg['hip']), fits.Column(name='tycho2_id', format='13A', array=sg['tycho2_id']), fits.Column(name='ra_tgas', format='D', array=sg['ra_tgas']), fits.Column(name='ra_error', format='D', array=sg['ra_error']), fits.Column(name='dec_tgas', format='D', array=sg['dec_tgas']), fits.Column(name='dec_error', format='D', array=sg['dec_error']), fits.Column(name='parallax', format='D', array=sg['parallax']), fits.Column(name='parallax_error', format='D', array=sg['parallax_error']), fits.Column(name='pmra', format='D', array=sg['pmra']), fits.Column(name='pmra_error', format='D', array=sg['pmra_error']), fits.Column(name='pmdec', format='D', array=sg['pmdec']), fits.Column(name='pmdec_error', format='D', array=sg['pmdec_error']), fits.Column(name='ra_dec_corr', format='E', array=sg['ra_dec_corr']), fits.Column(name='ra_parallax_corr', format='D', array=sg['ra_parallax_corr']), fits.Column(name='ra_pmra_corr', format='D', array=sg['ra_pmra_corr']), fits.Column(name='ra_pmdec_corr', format='D', array=sg['ra_pmdec_corr']), fits.Column(name='dec_parallax_corr', format='D', array=sg['dec_parallax_corr']), fits.Column(name='dec_pmra_corr', format='D', array=sg['dec_pmra_corr']), fits.Column(name='dec_pmdec_corr', format='D', array=sg['dec_pmdec_corr']), fits.Column(name='parallax_pmra_corr', format='D', array=sg['parallax_pmra_corr']), fits.Column(name='parallax_pmdec_corr', format='D', array=sg['parallax_pmdec_corr']), fits.Column(name='pmra_pmdec_corr', format='D', array=sg['pmra_pmdec_corr']), fits.Column(name='phot_g_mean_mag', format='D', array=sg['phot_g_mean_mag']), fits.Column(name='phot_g_mean_flux', format='D', array=sg['phot_g_mean_flux']), fits.Column(name='phot_g_mean_flux_error', format='D', array=sg['phot_g_mean_flux_error']), fits.Column(name='phot_variable_flag', format='19A', array=sg['phot_variable_flag']), fits.Column(name='gl_tgas', format='D', array=sg['gl_tgas']), fits.Column(name='gb_tgas', format='D', array=sg['gb_tgas']), fits.Column(name='ecl_lon', format='D', array=sg['ecl_lon']), fits.Column(name='ecl_lat', format='D', array=sg['ecl_lat']), fits.Column(name='angsep', format='D', array=sg['angsep']), fits.Column(name='dist', format='D', array=pc), fits.Column(name='distmod', format='D', array=distmod), fits.Column(name='Mg', format='D', array=Mg), fits.Column(name='ebv', format='D', array=ebv), fits.Column(name='parallax_hogg', format='D', array=parallax)])

    '''
    names = ['nuv_mag', 'fuv_mag', 'ra_gais', 'dec_gais', 'gl_gais', 'gb_gais', 'tilenum', 'Tycho_ID', 'Tycho_RA', 'Tycho_Dec', 'Tycho_Galactic_l', 'Tycho_Galactic_b', 'Tycho_B_Mag', 'Tycho_B_Mag_Err', 'Tycho_V_Mag', 'Tycho_V_Mag_Err', 'EBV', 'TMASS_ID', 'TMASS_J_Mag', 'TMASS_J_Mag_Err', 'TMASS_H_Mag', 'TMASS_H_Mag_Err', 'TMASS_K_Mag', 'TMASS_K_Mag_Err', 'Wise_ID', 'Wise_W1_Mag', 'Wise_W1_Mag_Err', 'Wise_W2_Mag', 'Wise_W2_Mag_Err', 'Wise_W3_Mag', 'Wise_W3_Mag_Err', 'Wise_W4_Mag', 'Wise_W4_Mag_Err', 'APASS_B_Mag', 'APASS_B_Mag_Err', 'APASS_V_Mag', 'APASS_V_Mag_Err', 'APASS_g_Mag', 'APASS_g_Mag_Err', 'APASS_r_Mag', 'APASS_r_Mag_Err', 'APASS_i_Mag', 'APASS_i_Mag_Err', 'APASS_B_V', 'APASS_B_V_Err', 'HIP_ID', 'HIP_Parallax', 'HIP_Parallax_Err', 'HIP_Proper_Motion_RA', 'HIP_Proper_Motion_RA_Err', 'HIP_Proper_Motion_Dec', 'HIP_Proper_Motion_Dec_Err', 'HIP_Hp_Mag', 'HIP_Hp_Mag_Err', 'HIP_B_V', 'HIP_V_I', 'APOGEE_Teff', 'APOGEE_Logg', 'APOGEE_FeH', 'APOGEE_AlphaFe', 'APOGEE_Heliocentric_RV', 'APOGEE_Heliocentric_RV_Err', 'RAVE_Teff', 'RAVE_Logg', 'RAVE_FeH', 'RAVE_AlphaFe', 'RAVE_Heliocentric_RV', 'RAVE_Heliocentric_RV_Err', 'RAVE_OFe', 'RAVE_MgFe', 'RAVE_AlFe', 'RAVE_SiFe', 'RAVE_CaFe', 'RAVE_NiFe', 'Pastel_Teff', 'Pastel_Logg', 'Pastel_FeH', 'LAMOST_Teff', 'LAMOST_Logg', 'LAMOST_FeH', 'LAMOST_AlphaFe', 'LAMOST_CFe', 'LAMOST_NFe', 'LAMOST_Mass', 'LAMOST_Age', 'LAMOST_Heliocentric_RV', 'LAMOST_Heliocentric_RV_Err', 'GALAH_Teff', 'GALAH_Logg', 'GALAH_FeH', 'GALAH_AlphaFe', 'GALAH_Heliocentric_RV', 'Average_Teff', 'Average_Logg', 'Average_FeH', 'Average_AlphaFe', 'Average_Heliocentric_RV', 'parallax', 'parallax_error', 'Gaia_Proper_Motion_RA', 'Gaia_Proper_Motion_RA_Err', 'Gaia_Proper_Motion_Dec', 'Gaia_Proper_Motion_Dec_Err', 'phot_g_mean_mag', 'Position_X', 'Position_Y', 'Position_Z', 'Velocity_X', 'Velocity_Y', 'Velocity_Z', 'Close_Pairs', 'angsep']

    formats = ['D', 'D', 'D', 'D', 'D', 'D', 'J', 'K', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'K', 'D', 'D', 'D', 'D', 'D', 'D', 'K', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'J', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'E', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'E', 'D']

    cols = []

    for column in range(len(names)):
        cols.append(fits.Column(name=names[column], format=formats[column], array=sg[names[column]]))
    cols.append(fits.Column(name='dist', format='D', array=pc))
    cols.append(fits.Column(name='distmod', format='D', array=distmod))
    cols.append(fits.Column(name='Mg', format='D', array=Mg))
    cols.append(fits.Column(name='ebv', format='D', array=ebv))
    cols.append(fits.Column(name='parallax_hogg', format='D', array=parallax))
    cols = fits.ColDefs(cols)
    '''

    endtable = fits.BinTableHDU.from_columns(cols)
    #endtable.writeto('../gais_tgasmatch_dust.fits')
    endtable.writeto('../gais_tgas_dust.fits')
