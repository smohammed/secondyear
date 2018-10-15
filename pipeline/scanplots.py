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


# For EXPmap
regions = ['0-10', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']

#incregions = ['0095', '0149', '0797', '0905', '0914', '1031', '1040', '1229', '1274', '2237', '2732', '2831', '2894', '3065', '3092', '3245', '3299', '3380', '3398', '3425', '3434', '3452', '3488', '3497', '3506', '3515', '3524', '3533', '3542', '3551', '3560', '3578']


for area in range(len(regions)):
    curregion = regions[area]
    print 'current region = ', curregion

    #exposure = Table()
    #expmap = fits.open('../../galexscans/count_map_'+curregion+'_exp.fits')[0].data
    
    counts = Table()
    ctmap = fits.open('../../galexscans/count_map_'+curregion+'_count.fits')[0].data
    #bkgd = fits.open('../Dunmaps/background/background_im1_'+curregion+'.fits')[0].data

    # Add WCS coordinates to get image counts
    hdulist = fits.open('../../galexscans/count_map_'+curregion+'_count.fits')
    w = wcs.WCS(hdulist[0].header)
    #pixels = np.array([[0, np.shape(expmap)[1]], [0, np.shape(expmap)[0]]]).T
    pixels = np.array([[0, np.shape(ctmap)[1]], [0, np.shape(ctmap)[0]]]).T
    world = w.wcs_pix2world(pixels, 1)

    gblimits = [world[0][1], world[1][1]]
    #gbrange = np.linspace(gblimits[0], gblimits[1], (np.shape(expmap)[0]))
    gbrange = np.linspace(gblimits[0], gblimits[1], (np.shape(ctmap)[0]))

    '''
    expvals = []
    for line in range(len(expmap[:, 0])):
        a = expmap[line,:][np.where(expmap[line,:] > 0.)]
        expvals.append(np.mean(a))

    # Add values to list
    exposure['gb'] = gbrange
    exposure[curregion] = expvals
    '''
    ctvals = []
    for line in range(len(ctmap[:, 0])):
        a = ctmap[line,:][np.where(ctmap[line,:] > 0.)]
        ctvals.append(np.mean(a))

    # Add values to list
    counts['gb'] = gbrange
    counts[curregion] = ctvals



    #ascii.write(exposure, '../../galexscans/exposuremaps/expmapvals'+curregion+'.txt', format='basic')
    ascii.write(counts, '../../galexscans/countmaps/ctmapvals'+curregion+'.txt', format='basic')
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

'''
for scan in regions:
    print scan
    a = Table.read('expmapvals'+scan+'.txt', format='ascii')
    a.rename_column(scan, 'expval')
    plt.plot(gbrange, a['expval'])

plt.xlabel('Galactic Latitude')
plt.ylabel('Average exposure time per scan')
plt.show()
'''

for scan in regions:
    print scan
    a = Table.read('ctmapvals'+scan+'.txt', format='ascii')
    a.rename_column(scan, 'ctval')
    plt.plot(gbrange, a['ctval'])

plt.xlabel('Galactic Latitude')
plt.ylabel('Average counts per scan')
plt.show()


