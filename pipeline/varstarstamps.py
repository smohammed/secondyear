from astropy.io import fits, ascii
import numpy as np
import os
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 17

regions = ['0.5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '110.3', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

#var = Table.read('../Dunmaps/fwhm/varstar_sex_gais.txt', format='ascii')
gv = Table.read('../Dunmaps/fwhm/varstar_mast_matches.txt', format='ascii')

gvgal = SkyCoord(gv['ra']*u.deg, gv['dec']*u.deg, frame='icrs').galactic

for i in range(len(gvgal)):
    for j in range(len(regions)):
        region = regions[j]
        if np.abs(gvgal.l.degree[i] - float(regions[j])) < 0.9:
            l, b = gvgal.l.degree[i], gvgal.b.degree[i]
            if not os.path.isfile('../images/varstamps/var_'+str(gv['uploadID'][i])+'_gl_'+str(l)[:5]+'.png'):
                if region == '0.5':
                    region = '5'
                    img = fits.open('../Dunmaps/countmaps/count_map_name_'+region+'_gal_sec_in.fits')
                else:
                    img = fits.open('../Dunmaps/countmaps/count_map_name_'+region.replace('.', '')+'_gal_sec_in.fits')

                w = wcs.WCS(img[0].header)
                xpix, ypix = w.all_world2pix(gvgal.l.degree[i], gvgal.b.degree[i], 1)

                offset = 100
                xmin = xpix - offset
                xmax = xpix + offset
                ymin = ypix - offset
                ymax = ypix + offset

                l, b = gvgal.l.degree[i], gvgal.b.degree[i]
                dx = w.all_pix2world(xpix-offset, ypix-offset, 1)[0] - w.all_pix2world(xpix, ypix, 1)[0]
                dy = w.all_pix2world(xpix, ypix, 1)[1] - w.all_pix2world(xpix-offset, ypix-offset, 1)[1]

                img = fits.open('../Dunmaps/countmaps/count_map_name_'+region.replace('.', '')+'_gal_sec_in.fits')[0].data
                img = img[ymin:ymax, xmin:xmax]
                plt.imshow(img, vmin=0, vmax=0.5, origin='lower', interpolation='nearest', aspect='auto', cmap=cm.gray, extent=[l-dx, l+dx, b-dy, b+dy])
                #plt.show()
                plt.savefig('../images/varstamps/var_'+str(gv['uploadID'][i])+'_gl_'+str(l)[:5]+'.png')
                plt.close()
                print 'j = ', j
    print 'i = ', i
