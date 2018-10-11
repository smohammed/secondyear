from astropy.io import fits, ascii
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


#regions = ['0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104', '0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194', '0203', '0212', '0221', '0230', '0239', '0248', '0257', '0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356', '0392', '0428', '0437', '0446', '0455', '0464', '0473', '0482', '0491', '0500', '0671', '0689', '0716', '0743', '0752', '0761', '0770', '0779', '0788', '0797', '0806', '0815', '0824', '0833', '0878', '0887', '0896', '0905', '0914', '0923', '0932', '0941', '0950', '0959', '0968', '0977', '0986', '0995', '1004', '1013', '1022', '1031', '1040', '1049', '1058', '1067', '1076', '1103', '1112', '1121', '1130', '1139', '1148', '1193', '1211', '1229', '1247', '1256', '1265', '1274', '1283', '1292', '1301', '1310', '1319', '1328', '1337', '1346', '1355', '1364', '1373', '1382', '1391', '1400', '1409', '1418', '1436', '1445', '1454', '1481', '1490', '1499', '1508', '1517', '1526', '1535', '1553', '1562', '1571', '1580', '1607', '1616', '1634', '1670', '1679', '1724', '1733', '1742', '1751', '1760', '1769', '1778', '1787', '1796', '1805', '1832', '1850', '1904', '1913', '1976', '1985', '2003', '2012', '2030', '2039', '2057', '2066', '2075', '2084', '2093', '2102', '2111', '2120', '2129', '2138', '2147', '2156', '2165', '2174', '2183', '2192', '2201', '2210', '2219', '2228', '2237', '2246', '2255', '2264', '2282', '2291', '2300', '2309', '2318', '2345', '2354', '2363', '2372', '2381', '2390', '2399', '2408', '2417', '2426', '2435', '2444', '2453', '2462', '2471', '2480', '2489', '2498', '2507', '2516', '2525', '2534', '2543', '2552', '2561', '2570', '2588', '2597', '2606', '2615', '2633', '2642', '2651', '2660', '2669', '2687', '2696', '2705', '2714', '2723', '2732', '2741', '2750', '2759', '2768', '2786', '2795', '2813', '2831', '2840', '2858', '2867', '2885', '2894', '2903', '2912', '2921', '2930', '2939', '2957', '2975', '2984', '3011', '3020', '3029', '3038', '3047', '3056', '3065', '3083', '3092', '3101', '3155', '3164', '3173', '3182', '3191', '3200', '3209', '3218', '3227', '3236', '3245', '3254', '3263', '3272', '3281', '3290', '3299', '3317', '3326', '3335', '3344', '3353', '3380', '3389', '3398', '3416', '3425', '3434', '3452', '3488', '3497', '3506', '3515', '3524', '3533', '3542', '3551', '3560', '3578', '3587', '3596']

regions = ['0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104', '0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194', '0203', '0212', '0221', '0230', '0239', '0248', '0257', '0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356', '0392', '0428', '0437', '0446', '0455', '0464', '0473', '0482', '0491', '0500', '0671', '0689', '0716', '0743', '0752', '0761', '0770', '0779', '0788', '0797', '0806', '0815', '0824', '0833', '0878', '0887', '0896', '0905', '0914', '0923', '0932', '0941', '0950', '0959', '0968', '0977', '0986', '0995', '1004', '1013', '1022', '1031', '1040', '1049', '1058', '1067', '1076', '1103', '1112', '1121', '1130', '1139', '1148', '1157', '1166', '1193', '1211', '1229', '1238', '1247', '1256', '1265', '1274', '1292', '1301', '1310', '1319', '1328', '1337', '1346', '1355', '1364', '1373', '1382', '1391', '1400', '1409', '1418', '1436', '1445', '1454', '1463', '1481', '1490', '1499', '1508', '1517', '1526', '1535', '1553', '1562', '1571', '1580', '1607', '1616', '1634', '1670', '1679', '1724', '1733', '1742', '1751', '1760', '1769', '1778', '1787', '1796', '1805', '1832', '1850', '1904', '1913', '1976', '1985', '2003', '2012', '1157', '1166', '1463', '2039', '2057', '2066', '2084', '2093', '2102', '2138', '2156', '2165', '2174', '2192', '2210', '2246', '2264', '2282', '2372', '2381', '2390', '2399', '2417', '2426', '2435', '2453', '2480', '2489', '2498', '2507', '2516', '2525', '2534', '2543', '2552', '2561', '2570', '2588', '2597', '2606', '2615', '2624', '2633', '2642', '2651', '2660', '2669', '2687', '2696', '2705', '2714', '2723', '2732', '2741', '2750', '2759', '2768', '2786', '2795', '2813', '2831', '2840', '2849', '2858', '2867', '2885', '2894', '2903', '2912', '2921', '2930', '2939', '2957', '2975', '2984', '3011', '3020', '3029', '3038', '3047', '3056', '3065', '3083', '3092', '3101', '3155', '3164', '3173', '3182', '3191', '3200', '3209', '3218', '3227', '3236', '3245', '3254', '3263', '3272', '3281', '3290', '3299', '3317', '3326', '3335', '3344', '3353', '3371', '3380', '3389', '3398', '3416', '3425', '3434', '3452', '3488', '3497', '3506', '3515', '3524', '3533', '3542', '3551', '3560', '3578', '3587', '3596', '1283', '1427', '1472', '1544','1589', '1598', '1652', '1661']

#incregions = ['0095', '0149', '0797', '0905', '0914', '1031', '1040', '1229', '1274', '2237', '2732', '2831', '2894', '3065', '3092', '3245', '3299', '3380', '3398', '3425', '3434', '3452', '3488', '3497', '3506', '3515', '3524', '3533', '3542', '3551', '3560', '3578']

#incregiondict = dict({'0095': [1214, 3950, 12318, 52037], '0149': [1214, 3950, 3532, 29730], '0797': [1214, 3950, 26493, 51563], '0905': [1214, 3950, 13600, 51230], '0914': [1214, 3950, 18800, 51600], '1031': [1214, 3950, 9300, 51600], '1040': [1214, 3950, 26400, 51600], '1229': [1214, 3950, 35000, 51600], '1274': [1214, 3950, 21800, 51300], '2237': [1214, 3950, 3300, 47000], '2732': [1214, 3950, 13100, 52300], '2831': [1214, 3950, 5700, 51600], '2894': [1214, 3950, 37000, 51500], '3065': [1214, 3950, 19000, 51500], '3092': [1214, 3950, 33400, 52400], '3245': [1214, 3950, 3200, 45000], '3299': [1214, 3950, 27600, 44500], '3380': [1214, 3950, 36000, 53000], '3398': [1605, 3490, 7200, 53000], '3425': [1214, 3950, 3000, 42900], '3434': [1214, 3950, 12200, 52600], '3452': [1214, 3950, 14600, 40800], '3488': [1214, 3950, 30400, 52600], '3497': [1214, 3950, 12100, 52600], '3506': [1214, 3950, 13200, 52600], '3515': [1214, 3950, 13380, 52600], '3524': [1214, 3950, 15500, 52700], '3533': [1214, 3950, 16400, 52700], '3542': [1214, 3950, 17000, 52700], '3551': [1214, 3950, 17700, 52700], '3560': [1214, 3950, 22700, 52700], '3578': [1214, 3950, 23600, 52700]})

#files = np.loadtxt('../scst/scst_list_firstscan.txt', dtype='str')


background = Table()
counts = Table()

for area in range(len(regions)):
    #scan = fits.open('../scst/'+files[area])[1].data
    curregion = regions[area]

    print 'current region = ', curregion

    exposure = Table()
    '''
    # Set correct image sizes
    if curregion in incregions:
        im1xmin, im1xmax, im1ymin, im1ymax = incregiondict[curregion][0], incregiondict[curregion][1], incregiondict[curregion][2], incregiondict[curregion][3]
    else:
        im1xmin, im1xmax, im1ymin, im1ymax = 1214, 3950, 3532, 51230
    '''
    # Load all tables and images, resize images
    #t1 = Table.read('../Dunmaps/fwhm/11-18data/starcat_'+curregion+'mapweight_fwhm.txt', format='ascii')
    #t1.remove_columns(('X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'x_new', 'y_new', 'FLUXERR_AUTO', 'FLUX_APER', 'ra', 'dec'))
    #t1 = t1[np.where((t1['FWHM_IMAGE'] < 10) & (t1['FWHM_IMAGE'] > 0))]

    #img = fits.open('../Dunmaps/countmaps/count_map_name_'+curregion+'_in.fits')[0].data
    #img = img[im1ymin:im1ymax, im1xmin:im1xmax]

    expmap = fits.open('../../galexscans/exposuremaps/count_map_'+curregion+'_exp.fits')[0].data
    #expmap = expmap[im1ymin:im1ymax, im1xmin:im1xmax]
    #bkgd = fits.open('../Dunmaps/background/background_im1_'+curregion+'.fits')[0].data

    # Add WCS coordinates to get image counts
    hdulist = fits.open('../../galexscans/exposuremaps/count_map_'+curregion+'_exp.fits')
    w = wcs.WCS(hdulist[0].header)
    #pixels = np.array([[im1xmin, im1xmax], [im1ymin, im1ymax]]).T
    pixels = np.array([[0, np.shape(expmap)[1]], [0, np.shape(expmap)[0]]]).T
    world = w.wcs_pix2world(pixels, 1)

    gblimits = [world[0][1], world[1][1]]
    gbrange = np.linspace(gblimits[0], gblimits[1], (np.shape(expmap)[0]))
    #gbrange = np.linspace(gblimits[0], gblimits[1], (im1ymax-im1ymin))

    # Add values to list
    expmapvals = expmap[:, np.shape(expmap)[1]/2]

    exposure['gb'] = gbrange
    exposure[curregion] = expmapvals

    ascii.write(exposure, '../../galexscans/exposuremaps/expmapvals'+curregion+'.txt', format='basic')
    '''
    bkgdvals = []
    for line in range(len(bkgd[:, 0])):
        bkgdvals.append(np.mean(bkgd[line, :]))
    bkgdvals = np.array(bkgdvals)

    imgvals = []
    for line in range(0, len(img[:, 0]), 20):
        imgvals.append(np.mean(img[line, :]))
    imgvals = np.array(imgvals)

    imggbrange = np.linspace(gblimits[0], gblimits[1], (im1ymax-im1ymin)/20)

    #if curregion in incregions:
    print len(imggbrange)
    print len(imgvals)

    if len(imggbrange) != len(imgvals):
        if len(imggbrange) > len(imgvals):
                imggbrange = np.delete(imggbrange, -1)
        elif len(imggbrange) < len(imgvals):
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


    exposure[curregion] = expmapvals
    background[curregion] = bkgdvals
    counts[curregion] = imgvals


    # Get histogram for number of sources
    n, bins = np.histogram(t1['gb'], bins=100)
    bins_mean = [0.5*(bins[i] + bins[i+1]) for i in range(len(n))]


    # Now plot
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    scangal = SkyCoord(scan['ra_acs']*u.deg, scan['dec_acs']*u.deg, frame='icrs').galactic
    ax1.plot(scangal.b.degree, scan['NDCTEC'], label='TEC')
    ax1.plot(scangal.b.degree, scan['NDCFEC'], label='FEC', c='red')

    ax1.set_xlim((-10, 10))
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0., fontsize=15)

    # exp, background and counts vs gl

    ax2.plot(imggbrange, imgvals*1000.-20, label='counts*1000', c='green')
    ax2.plot(avgimggbrange, avgimgvals*1000*2, label='avgcount*1000/deg', c='orange', linewidth=2)
    ax2.plot(gbrange, expmapvals, label='expmap', c='blue', linewidth=2)
    ax2.plot(gbrange, bkgdvals*4000. - 70., label='bkgd*4000', c='red', linewidth=2)
    ax2.plot(bins_mean, n/100.-20, label='# sources/100', c='black', linewidth=2)
    ax2.set_xlim((-10, 10))
    if np.max(expmapvals) > 120:
        ax2.set_ylim(-50, np.max(expmapvals) + 2)
    else:
        ax2.set_ylim(-50, 120)
    ax2.set_xlabel('gb, for scan '+files[area][13:23])
    ax2.set_ylabel('counts')
    ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0., fontsize=15)

    #plt.suptitle('region = '+curregion)
    plt.tight_layout()
    plt.savefig(curregion+'_scanplots.png')
    #plt.show()
    plt.close()
    '''

#ascii.write(background, 'bkgdmapvals.txt', format='basic')
#ascii.write(counts, 'countmapvals.txt', format='basic')








