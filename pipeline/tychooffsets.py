from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 17


# Input parameters
full = 1
partial = 0

fec = 1

#################################################
# Subplots for each scan
#################################################
# All complete scans
if (fec == 0) and (full == 1):
    scans = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '10.4', '11.3', '12.2', '14.0', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '124.7', '125.6', '126.5', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '284.0', '285.8', '286.7', '288.5', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '308.3', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '325.4', '326.3', '327.2', '328.1', '329.0', '331.7', '332.6', '333.5', '334.4', '335.3', '338.9', '341.6', '358.7', '359.6']

# fec scans
if (fec == 1) and (full == 1):
    #scans = ['0014', '0059', '0203', '0239', '1310', '1607', '3209', '3236']
    scans = ['0014', '0032', '0059', '0743']

# Incomplete scans
if (fec == 0) and (partial == 1):
    incscans = ['9.5', '14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']

#################################################
# Load catalogs
#################################################
tycho = Table(fits.open('../tycho2.fits')[1].data)
tycho.remove_columns(('RAJ2000', 'DEJ2000', 'TYC1', 'TYC2', 'TYC3', 'pmRA', 'pmDE', 'BTmag', 'e_BTmag', 'VTmag', 'e_VTmag', 'HIP', 'RA_ICRS_', 'DE_ICRS_'))
tychogal = SkyCoord(tycho['Glon']*u.deg, tycho['Glat']*u.deg, frame='galactic')

# Label partial or full scans
if full == 1:
    skyrange = scans
if partial == 1:
    skyrange = incscans

for curregion in skyrange:
    print curregion

    if fec == 0:
        t1 = Table.read('../Dunmaps/fwhm/03-02data/starcat_'+curregion+'mapweight_fwhm.txt', format='ascii')

    if fec == 1:
        t1 = Table.read('../Dunmaps/fwhm/fec/03-02data/starcat_'+curregion+'mapweight_fec_fwhm.txt', format='ascii')

    t1.remove_columns(('X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'x_new', 'y_new', 'FLUXERR_AUTO', 'FLUX_APER', 'ra', 'dec'))
    t1 = t1[np.where((t1['FWHM_IMAGE'] < 10) & (t1['FWHM_IMAGE'] > 3.5))]
    t1gal = SkyCoord(t1['gl']*u.deg, t1['gb']*u.deg, frame='galactic')

    # Match to Tycho 2
    t1tyind, tychoind, angsepty, ang3d = search_around_sky(t1gal, tychogal, 3*u.arcsec)
    t1ty = hstack([t1[t1tyind], Table(tycho[tychoind])])
    t1ty['angsep'] = angsepty
    t1ty.rename_column('gl', 'gl_sex')
    t1ty.rename_column('gb', 'gb_sex')
    t1ty.rename_column('Glon', 'gl_tycho')
    t1ty.rename_column('Glat', 'gb_tycho')

    if curregion == '5':
        t1ty['gl_sex'][np.where(t1ty['gl_sex'] > 350)] = t1ty['gl_sex'][np.where(t1ty['gl_sex'] > 350)] - 360
        t1ty['gl_tycho'][np.where(t1ty['gl_tycho'] > 350)] = t1ty['gl_tycho'][np.where(t1ty['gl_tycho'] > 350)] - 360

    if curregion == '359.6':
        t1ty['gl_sex'][np.where(t1ty['gl_sex'] < 350)] = t1ty['gl_sex'][np.where(t1ty['gl_sex'] < 350)] + 360
        t1ty['gl_tycho'][np.where(t1ty['gl_tycho'] < 350)] = t1ty['gl_tycho'][np.where(t1ty['gl_tycho'] < 350)] + 360

    print 'SExtractor catalog loaded, matched with Tycho2 and GAIS'

    #################################################
    # Plot everything
    #################################################
    fig = plt.figure()
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((3, 2), (1, 0), colspan=2)
    ax3 = plt.subplot2grid((3, 2), (2, 0))
    ax4 = plt.subplot2grid((3, 2), (2, 1))

    # T2-sex gl
    sc1 = ax1.scatter(t1ty['gl_tycho'], (t1ty['gl_sex']-t1ty['gl_tycho'])*3600, edgecolor='none', c=t1ty['gb_tycho'], vmin=-10, vmax=10)
    ax1.set_xlim((np.min(t1ty['gl_tycho']), np.max(t1ty['gl_tycho'])))
    ax1.set_ylim((-3, 3))
    ax1.set_xlabel('gl (T2) [deg]')
    ax1.set_ylabel('gl (SEx-T2) [arcsec]')
    cbar1 = plt.colorbar(sc1, ax=ax1)
    cbar1.ax.tick_params(labelsize=12)
    cbar1.set_label('gb (T2)', size=14)

    # T2-sex gb
    sc2 = ax2.scatter(t1ty['gl_tycho'], (t1ty['gb_sex']-t1ty['gb_tycho'])*3600, edgecolor='none', c=t1ty['gb_tycho'], vmin=-10, vmax=10)  # vmin=np.min(t1ty['gb_tycho']), vmax=np.max(t1ty['gb_tycho']))
    ax2.set_xlim((np.min(t1ty['gl_tycho']), np.max(t1ty['gl_tycho'])))
    ax2.set_ylim((-3, 3))
    ax2.set_xlabel('gl (T2) [deg]')
    ax2.set_ylabel('gb (SEx-T2) [arcsec]')
    cbar2 = plt.colorbar(sc2, ax=ax2)
    cbar2.ax.tick_params(labelsize=12)
    cbar2.set_label('gb (T2)', size=14)

    # Tycho match distance
    ax3.scatter((t1ty['gl_sex']-t1ty['gl_tycho'])*3600, (t1ty['gb_sex']-t1ty['gb_tycho'])*3600, edgecolor='none', alpha=0.2)
    ax3.axhline(y=0, c='red')
    ax3.axvline(x=0, c='red')
    ax3.set_xlim((-3, 3))
    ax3.set_ylim((-3, 3))
    ax3.set_xlabel('gl (SEx-T2), N = '+str(len(t1ty)))
    ax3.set_ylabel('gb (SEx-T2)')

    # NUV histogram
    area = 20 * 0.9  # 6118 # 360*20 - blank area
    x, bins, p = ax4.hist(t1ty['nuv'], bins=np.linspace(12, 25, 25))
    for item in p:
        item.set_height(np.log10(item.get_height()/area))
    ax4.set_xlim((13, 20))
    ax4.set_ylim(0, 1.5)
    ax4.set_xlabel('NUV (T2 + SEx)')
    ax4.set_ylabel('counts/area')
    if fec == 0:
        if curregion == '5':
            plt.suptitle('partial scan, region = 0.5, len = '+str(len(t1ty)))
        else:
            plt.suptitle('partial scan, region = '+curregion+', len = '+str(len(t1ty)))
    if fec == 1:
        plt.suptitle('FEC region = '+curregion+', len = '+str(len(t1ty)))

    plt.tight_layout()

    if fec == 0:
        if curregion == '5':
            plt.savefig('../11-22-region0.5tychoplots.png')
        else:
            plt.savefig('../11-22-region'+curregion+'tychoplots.png')
            #images/11-22-tychoplots/
    if fec == 1:
        plt.savefig('../03-02-region'+curregion+'tychoplots_fec.png')

    #plt.show()
    plt.close()
