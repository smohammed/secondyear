from astropy.io import fits
from astropy.table import Table, hstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky


def catmatch(catalog, match):
    cat = fits.open('../'+catalog)[1].data
    catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')

    if match == 'gaia':
        gaia = fits.open('../gaia-tycho.fits')[1].data
        gaiagal = SkyCoord(gaia['l']*u.deg, gaia['b']*u.deg, frame='galactic')

        catind, gind, angsep, ang3d = search_around_sky(catgal, gaiagal, 3*u.arcsec)
        comb = hstack([Table(cat[catind]), Table(gaia[gind])])
        comb.rename_column('nuv', 'nuv_mag')
        comb.rename_column('gl', 'gl_sex')
        comb.rename_column('gb', 'gb_sex')
        comb.rename_column('ra', 'ra_sex')
        comb.rename_column('dec_1', 'dec_sex')
        comb.rename_column('ra_2', 'ra_gaia')
        comb.rename_column('dec_2', 'dec_gaia')
        comb.rename_column('l', 'gl_gaia')
        comb.rename_column('b', 'gb_gaia')
        comb['angsep_sg'] = angsep

        cols = fits.ColDefs([fits.Column(name='X_IMAGE', format='D', array=comb['X_IMAGE']), fits.Column(name='Y_IMAGE', format='D', array=comb['Y_IMAGE']), fits.Column(name='FLUX_APER', format='D', array=comb['FLUX_APER']), fits.Column(name='A_IMAGE', format='D', array=comb['A_IMAGE']), fits.Column(name='B_IMAGE', format='D', array=comb['B_IMAGE']), fits.Column(name='THETA_IMAGE', format='D', array=comb['THETA_IMAGE']), fits.Column(name='FWHM_IMAGE', format='D', array=comb['FWHM_IMAGE']), fits.Column(name='nuv_mag', format='D', array=comb['nuv_mag']), fits.Column(name='gl_sex', format='D', array=comb['gl_sex']), fits.Column(name='gb_sex', format='D', array=comb['gb_sex']), fits.Column(name='ra_sex', format='D', array=comb['ra_sex']), fits.Column(name='dec_sex', format='D', array=comb['dec_sex']), fits.Column(name='hip', format='K', array=comb['hip']), fits.Column(name='tycho2_id', format='13A', array=comb['tycho2_id']), fits.Column(name='solution_id', format='K', array=comb['solution_id']), fits.Column(name='source_id', format='K', array=comb['source_id']), fits.Column(name='random_index', format='K', array=comb['random_index']), fits.Column(name='ref_epoch', format='D', array=comb['ref_epoch']), fits.Column(name='ra_gaia', format='D', array=comb['ra_gaia']), fits.Column(name='ra_error', format='D', array=comb['ra_error']), fits.Column(name='dec_gaia', format='D', array=comb['dec_gaia']), fits.Column(name='dec_error', format='D', array=comb['dec_error']), fits.Column(name='parallax', format='D', array=comb['parallax']), fits.Column(name='parallax_error', format='D', array=comb['parallax_error']), fits.Column(name='pmra', format='D', array=comb['pmra']), fits.Column(name='pmra_error', format='D', array=comb['pmra_error']), fits.Column(name='pmdec', format='D', array=comb['pmdec']), fits.Column(name='pmdec_error', format='D', array=comb['pmdec_error']), fits.Column(name='ra_dec_corr', format='E', array=comb['ra_dec_corr']), fits.Column(name='ra_parallax_corr', format='D', array=comb['ra_parallax_corr']), fits.Column(name='ra_pmra_corr', format='D', array=comb['ra_pmra_corr']), fits.Column(name='ra_pmdec_corr', format='D', array=comb['ra_pmdec_corr']), fits.Column(name='dec_parallax_corr', format='D', array=comb['dec_parallax_corr']), fits.Column(name='dec_pmra_corr', format='D', array=comb['dec_pmra_corr']), fits.Column(name='dec_pmdec_corr', format='D', array=comb['dec_pmdec_corr']), fits.Column(name='parallax_pmra_corr', format='D', array=comb['parallax_pmra_corr']), fits.Column(name='parallax_pmdec_corr', format='D', array=comb['parallax_pmdec_corr']), fits.Column(name='pmra_pmdec_corr', format='D', array=comb['pmra_pmdec_corr']), fits.Column(name='astrometric_n_obs_al', format='K', array=comb['astrometric_n_obs_al']), fits.Column(name='astrometric_n_obs_ac', format='K', array=comb['astrometric_n_obs_ac']), fits.Column(name='astrometric_n_good_obs_al', format='K', array=comb['astrometric_n_good_obs_al']), fits.Column(name='astrometric_n_good_obs_ac', format='K', array=comb['astrometric_n_good_obs_ac']), fits.Column(name='astrometric_n_bad_obs_al', format='K', array=comb['astrometric_n_bad_obs_al']), fits.Column(name='astrometric_n_bad_obs_ac', format='K', array=comb['astrometric_n_bad_obs_ac']), fits.Column(name='astrometric_delta_q', format='D', array=comb['astrometric_delta_q']), fits.Column(name='astrometric_excess_noise', format='D', array=comb['astrometric_excess_noise']), fits.Column(name='astrometric_excess_noise_sig', format='D', array=comb['astrometric_excess_noise']), fits.Column(name='astrometric_primary_flag', format='25A', array=comb['astrometric_primary_flag']), fits.Column(name='astrometric_relegation_factor', format='D', array=comb['astrometric_relegation_factor']), fits.Column(name='astrometric_weight_al', format='D', array=comb['astrometric_weight_al']), fits.Column(name='astrometric_weight_ac', format='D', array=comb['astrometric_weight_ac']), fits.Column(name='astrometric_priors_used', format='K', array=comb['astrometric_priors_used']), fits.Column(name='matched_observations', format='K', array=comb['matched_observations']), fits.Column(name='duplicated_source', format='18A', array=comb['duplicated_source']), fits.Column(name='scan_direction_strength_k1', format='D', array=comb['scan_direction_strength_k1']), fits.Column(name='scan_direction_strength_k2', format='D', array=comb['scan_direction_strength_k2']), fits.Column(name='scan_direction_strength_k3', format='D', array=comb['scan_direction_strength_k3']), fits.Column(name='scan_direction_strength_k4', format='D', array=comb['scan_direction_strength_k4']), fits.Column(name='scan_direction_mean_k1', format='D', array=comb['scan_direction_mean_k1']), fits.Column(name='scan_direction_mean_k2', format='D', array=comb['scan_direction_mean_k2']), fits.Column(name='scan_direction_mean_k3', format='D', array=comb['scan_direction_mean_k3']), fits.Column(name='scan_direction_mean_k4', format='D', array=comb['scan_direction_mean_k4']), fits.Column(name='phot_g_n_obs', format='J', array=comb['phot_g_n_obs']), fits.Column(name='phot_g_mean_flux', format='D', array=comb['phot_g_mean_flux']), fits.Column(name='phot_g_mean_flux_error', format='D', array=comb['phot_g_mean_flux_error']), fits.Column(name='phot_g_mean_mag', format='D', array=comb['phot_g_mean_mag']), fits.Column(name='phot_variable_flag', format='19A', array=comb['phot_variable_flag']), fits.Column(name='gl_gaia', format='D', array=comb['gl_gaia']), fits.Column(name='gb_gaia', format='D', array=comb['gb_gaia']), fits.Column(name='ecl_lon', format='D', array=comb['ecl_lon']), fits.Column(name='ecl_lat', format='D', array=comb['ecl_lat']), fits.Column(name='angsep_sg', format='D', array=comb['angsep_sg'])])

        endtable = fits.BinTableHDU.from_columns(cols)
        endtable.writeto('../sex_gaia-09-14.fits')

    if match == 'vphas':
        vphas = fits.open('../vphas_allg.fits')[1].data
        vpgal = SkyCoord(vphas['RAJ2000']*u.deg, vphas['DEJ2000']*u.deg, frame='icrs').galactic

        catind, vpind, angsep, ang3d = search_around_sky(catgal, vpgal, 3*u.arcsec)
        comb = hstack([Table(cat[catind]), Table(vphas[vpind])])
        comb.rename_column('nuv', 'nuv_mag')
        comb.rename_column('gl', 'gl_sex')
        comb.rename_column('gb', 'gb_sex')
        comb.rename_column('ra', 'ra_sex')
        comb.rename_column('dec', 'dec_sex')
        comb.rename_column('RAJ2000', 'ra_vphas')
        comb.rename_column('DEJ2000', 'dec_vphas')
        comb['gl_vphas'] = vpgal.l.degree[vpind]
        comb['gb_vphas'] = vpgal.b.degree[vpind]
        comb['angsep_sv'] = angsep

        cols = fits.ColDefs([fits.Column(name='X_IMAGE', format='D', array=comb['X_IMAGE']), fits.Column(name='Y_IMAGE', format='D', array=comb['Y_IMAGE']), fits.Column(name='FLUX_APER', format='D', array=comb['FLUX_APER']), fits.Column(name='A_IMAGE', format='D', array=comb['A_IMAGE']), fits.Column(name='B_IMAGE', format='D', array=comb['B_IMAGE']), fits.Column(name='THETA_IMAGE', format='D', array=comb['THETA_IMAGE']), fits.Column(name='FWHM_IMAGE', format='D', array=comb['FWHM_IMAGE']), fits.Column(name='nuv_mag', format='D', array=comb['nuv_mag']), fits.Column(name='gl_sex', format='D', array=comb['gl_sex']), fits.Column(name='gb_sex', format='D', array=comb['gb_sex']), fits.Column(name='ra_sex', format='D', array=comb['ra_sex']), fits.Column(name='dec_sex', format='D', array=comb['dec_sex']), fits.Column(name='u_AB', format='D', array=comb['u_AB']), fits.Column(name='g_AB', format='D', array=comb['g_AB']), fits.Column(name='r_AB', format='D', array=comb['r_AB']), fits.Column(name='r2_AB', format='D', array=comb['r2_AB']), fits.Column(name='i_AB', format='D', array=comb['i_AB']), fits.Column(name='angsep_sv', format='D', array=comb['angsep_sv'])])

        endtable = fits.BinTableHDU.from_columns(cols)
        endtable.writeto('../sex_vphas.fits')


catmatch('starcat_fwhm_12-14.fits', 'vphas')
