from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

# Input parameters
full = 1
partial = 0

fec = 1

#################################################
# Subplots for each scan
#################################################
# All complete scans
#scans = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '10.4', '11.3', '12.2', '14.0', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '124.7', '125.6', '126.5', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5'

#scans = ['253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '284.0', '285.8', '286.7', '288.5', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '308.3', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '325.4', '326.3', '327.2', '328.1', '329.0', '331.7', '332.6', '333.5', '334.4', '335.3', '338.9', '341.6', '358.7', '359.6']

# fec regions 04/17
#scans = ['0014', '0032', '0059', '0203', '0239']

# regions from 07/19
scans = ['0023', '0239', '0032', '0203', '0446', '0464', '0473', '0806', '0815', '1301', '1310', '1319', '1616', '1634', '1679', '2174', '2183', '2192', '2714', '2750', '3236', '3245', '3281']

# Incomplete scans
incscans = ['9.5', '14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']
incscandict = dict({'9.5': [1214, 3950, 12318, 52037], '14.9': [1214, 3950, 3532, 29730], '79.7': [1214, 3950, 26493, 51563], '90.5': [1214, 3950, 13600, 51230], '91.4': [1214, 3950, 18800, 51600], '103.1': [1214, 3950, 9300, 51600], '104.0': [1214, 3950, 26400, 51600], '122.9': [1214, 3950, 35000, 51600], '127.4': [1214, 3950, 21800, 51300], '223.7': [1214, 3950, 3300, 47000], '273.2': [1214, 3950, 13100, 52300], '283.1': [1214, 3950, 5700, 51600], '289.4': [1214, 3950, 37000, 51500], '306.5': [1214, 3950, 19000, 51500], '309.2': [1214, 3950, 33400, 52400], '324.5': [1214, 3950, 3200, 45000], '329.9': [1214, 3950, 27600, 44500], '338.0': [1214, 3950, 36000, 53000], '339.8': [1605, 3490, 7200, 53000], '342.5': [1214, 3950, 3000, 42900], '343.4': [1214, 3950, 12200, 52600], '345.2': [1214, 3950, 14600, 40800], '348.8': [1214, 3950, 30400, 52600], '349.7': [1214, 3950, 12100, 52600], '350.6': [1214, 3950, 13200, 52600], '351.5': [1214, 3950, 13380, 52600], '352.4': [1214, 3950, 15500, 52700], '353.3': [1214, 3950, 16400, 52700], '354.2': [1214, 3950, 17000, 52700], '355.1': [1214, 3950, 17700, 52700], '356.0': [1214, 3950, 22700, 52700], '357.8': [1214, 3950, 23600, 52700]})

#################################################
# Load catalogs
#################################################
tycho = Table(fits.open('../tycho2.fits')[1].data)
tycho.remove_columns(('RAJ2000', 'DEJ2000', 'TYC1', 'TYC2', 'TYC3', 'pmRA', 'pmDE', 'BTmag', 'e_BTmag', 'VTmag', 'e_VTmag', 'HIP', 'RA_ICRS_', 'DE_ICRS_'))
galex = Table(fits.open('../GALEXPlane2_smohammed.fit')[1].data)
galex.remove_columns(('fuv_mag', 'ra', 'dec'))
galexcut = np.where(galex['nuv_mag'] == -999.)
galex.remove_rows(galexcut)
tychogal = SkyCoord(tycho['Glon']*u.deg, tycho['Glat']*u.deg, frame='galactic')
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')

# Label partial or full scans
if full == 1:
    skyrange = scans
if partial == 1:
    skyrange = incscans

for curregion in skyrange:
    if full == 1:
        im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 3532, 51230
    if partial == 1:
        im1xmin, im1xmax, im1ymin, im1ymax = incscandict[curregion][0], incscandict[curregion][1], incscandict[curregion][2], incscandict[curregion][3]

    print curregion

    t1 = Table.read('../Dunmaps/fwhm/fec/07-19data/starcat_'+curregion+'mapweight_fec_fwhm.txt', format='ascii')

    #t1 = t1[np.where((t1['FWHM_IMAGE'] < 10) & (t1['FWHM_IMAGE'] > 3.5) & (t1['FLUX_AUTO'] > 0))]
    t1 = t1[np.where((t1['FLUX_AUTO'] > 0))]

    t1gal = SkyCoord(t1['gl']*u.deg, t1['gb']*u.deg, frame='galactic')
    t1tyind, tychoind, angsepty, ang3d = search_around_sky(t1gal, tychogal, 3*u.arcsec)
    t1ty = hstack([t1[t1tyind], Table(tycho[tychoind])])
    t1ty['angsep'] = angsepty
    t1ty.rename_column('gl', 'gl_sex')
    t1ty.rename_column('gb', 'gb_sex')
    t1ty.rename_column('Glon', 'gl_tycho')
    t1ty.rename_column('Glat', 'gb_tycho')

    t1gaisind, galexind, angsepga, ang3d = search_around_sky(t1gal, galexgal, 3*u.arcsec)
    t1gais = hstack([t1[t1gaisind], Table(galex[galexind])])
    t1gais['angsep'] = angsepga
    t1gais.rename_column('nuv', 'nuv_sex')
    t1gais.rename_column('nuv_mag', 'nuv_galex')

    print 'SExtractor catalog loaded, matched with Tycho2 and GAIS'

    '''
    #################################################
    # Load maps
    #################################################
    hdulist = fits.open('../fecmaps/07-19/count_map_'+curregion.replace('.', '')+'-cal-sec_in_dis_new_bp.fits')
    img = hdulist[0].data
    expmap = fits.open('../fecmaps/07-19/count_map_'+curregion.replace('.', '')+'-cal-sec_exp_bp.fits')[0].data
    bkgd = fits.open('../Dunmaps/background_im1_'+curregion.replace('.', '')+'_fec.fits')[0].data

    img = img[im1ymin:im1ymax, im1xmin:im1xmax]
    expmap = expmap[im1ymin:im1ymax, im1xmin:im1xmax]

    w = wcs.WCS(hdulist[0].header)
    pixels = np.array([[im1xmin, im1xmax], [im1ymin, im1ymax]]).T
    world = w.wcs_pix2world(pixels, 1)

    gblimits = [world[0][1], world[1][1]]
    gbrange = np.linspace(gblimits[0], gblimits[1], (im1ymax-im1ymin))

    expmapvals = expmap[:, np.shape(expmap)[1]/2]

    bkgdvals = []
    for line in range(len(bkgd[:, 0])):
        bkgdvals.append(np.mean(bkgd[line, :]))
    bkgdvals = np.array(bkgdvals)

    imgvals = []
    for line in range(0, len(img[:, 0]), 20):
        imgvals.append(np.mean(img[line, :]))
    imgvals = np.array(imgvals)
    imggbrange = np.linspace(gblimits[0], gblimits[1], (im1ymax-im1ymin)/20)
    if len(imgvals) != len(imggbrange):
        imgvals = np.delete(imgvals, -1)

    imgsize = (world[1][1] - world[0][1]) / len(img)  # About 19.6deg/pix, or 1.48arcsec/pix
    pixperdeg = 1/imgsize/2  # about 2400 pixels per degree

    # Make avg value for img over 0.5 degrees
    pixrange = []
    for line in np.arange(0, len(img[:, 0]), pixperdeg):
        pixrange.append(line)

    avgimgvals = []
    for i in range(len(pixrange)-1):
        avgimgvals.append(np.average(img[pixrange[i]:pixrange[i+1], :]))

    avgimgvals.append(np.average(img[pixrange[-1]:len(img), :]))
    avgimgvals = np.array(avgimgvals)

    avgimggbrange = np.linspace(gblimits[0], gblimits[1], (len(pixrange)-1)) - 0.5
    avgimggbrange = np.append(avgimggbrange, avgimggbrange[-1]+0.5)

    if curregion == '5':
        for i in range(len(t1ty)):
            if t1ty['gl_tycho'][i] > 350:
                t1ty['gl_tycho'][i] = t1ty['gl_tycho'][i] - 360

    print 'maps loaded'
    '''

    #################################################
    # Make histogram for sources
    #################################################
    n, bins = np.histogram(t1['gb'], bins=100)
    bins_mean = [0.5*(bins[i] + bins[i+1]) for i in range(len(n))]

    #################################################
    # Plot everything
    #################################################
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (0, 1))
    ax3 = plt.subplot2grid((2, 2), (1, 0))
    ax4 = plt.subplot2grid((2, 2), (1, 1))
    #ax5 = plt.subplot2grid((2, 2), (2, 0), colspan=2)

    # Tycho match distance
    ax1.scatter((t1ty['gl_sex']-t1ty['gl_tycho'])*3600, (t1ty['gb_sex']-t1ty['gb_tycho'])*3600, edgecolor='none', alpha=0.2)
    ax1.axhline(y=0, c='red')
    ax1.axvline(x=0, c='red')
    ax1.set_xlim((-3, 3))
    ax1.set_ylim((-3, 3))
    ax1.set_title('gl (SEx-T2), N='+str(len(t1ty)))
    ax1.set_ylabel('gb (SEx-T2)')

    # GAIS vs SExtractor NUV comparison
    im = ax2.scatter(t1gais['nuv_galex'], t1gais['nuv_sex']-t1gais['nuv_galex'], alpha=0.5, edgecolor='none', s=1, c=t1gais['gb'], vmin=-10, vmax=10)
    try:
        cbar = plt.colorbar(im, ax=ax2)
        cbar.ax.tick_params(labelsize=12)
    except TypeError:
        print 'len t1gais = ', len(t1gais)
    ax2.axhline(y=0, c='red')
    ax2.set_xlim((13, 25))
    ax2.set_ylim((-1, 1.5))
    ax2.set_ylabel('NUV (GAIS)')
    ax2.set_ylabel('NUV (SEx-GAIS)')
    ax2.set_title('NUV (GAIS), N = '+str(len(t1gais)))

    # NUV histogram
    area = 20 * 0.9  # 6118 # 360*20 - blank area
    x, bins, p = ax3.hist(t1['nuv'], bins=np.linspace(12, 25, 25))
    for item in p:
        item.set_height(np.log10(item.get_height()/area))
    ax3.set_xlim((13, 25))
    ax3.set_ylim(0, 4)
    ax3.set_title('NUV (only SExtractor)')
    ax3.set_ylabel('counts/area')

    # FWHM vs NUV
    #scatter_contour(t1['nuv'], t1['FWHM_IMAGE'], threshold=700, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
    ax4.scatter(t1['nuv'], t1['FWHM_IMAGE'], edgecolor='none', c='k', s=1)
    ax4.axhline(y=0,  c='red')
    ax4.set_xlim((13, 25))
    ax4.set_ylim((4,  14))
    ax4.set_title('NUV (Only SExtractor),  N = '+str(len(t1)))
    ax4.set_ylabel('FWHM')

    '''
    print len(imggbrange)
    print len(imgvals)

    # exp, background and counts vs gl
    ax5.plot(imggbrange, imgvals*1000.-20, label='counts*1000', c='green')
    ax5.plot(avgimggbrange, avgimgvals*1000, label='avgcount*1000', c='orange', linewidth=2)
    ax5.plot(gbrange, expmapvals, label='expmap', c='blue', linewidth=2)
    ax5.plot(gbrange, bkgdvals*4000. - 70., label='bkgd*4000', c='red', linewidth=2)
    ax5.plot(bins_mean, n/100.-20, label='# sources/100', c='black', linewidth=2)
    ax5.set_xlim((-10, 10))
    if np.max(expmapvals) > 120:
        ax5.set_ylim(-50, np.max(expmapvals) + 2)
    else:
        ax5.set_ylim(-50, 120)
    ax5.set_xlabel('gb')
    ax5.set_ylabel('counts')
    ax5.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0., fontsize=15)
    '''

    if curregion == '5':
        plt.suptitle('region = 0.5')
    else:
        plt.suptitle('pixfix, '+curregion)
    plt.tight_layout()

    plt.savefig('../images/07-19-region'+curregion+'matchplots_fec_doublepix.png')

    #plt.show()
    plt.close()
