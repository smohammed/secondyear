from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import os
from astropy.convolution import convolve, Gaussian2DKernel
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20


regions = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '110.3', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']


incregions = ['9.5', '14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']
incregiondict = dict({'9.5': [1214, 3950, 12318, 52037], '14.9': [1214, 3950, 3532, 29730], '79.7': [1214, 3950, 26493, 51563], '90.5': [1214, 3950, 13600, 51230], '91.4': [1214, 3950, 18800, 51600], '103.1': [1214, 3950, 9300, 51600], '104.0': [1214, 3950, 26400, 51600], '122.9': [1214, 3950, 35000, 51600], '127.4': [1214, 3950, 21800, 51300], '223.7': [1214, 3950, 3300, 47000], '273.2': [1214, 3950, 13100, 52300], '283.1': [1214, 3950, 5700, 51600], '289.4': [1214, 3950, 37000, 51500], '306.5': [1214, 3950, 19000, 51500], '309.2': [1214, 3950, 33400, 52400], '324.5': [1214, 3950, 3200, 45000], '329.9': [1214, 3950, 27600, 44500], '338.0': [1214, 3950, 36000, 53000], '339.8': [1605, 3490, 7200, 53000], '342.5': [1214, 3950, 3000, 42900], '343.4': [1214, 3950, 12200, 52600], '345.2': [1214, 3950, 14600, 40800], '348.8': [1214, 3950, 30400, 52600], '349.7': [1214, 3950, 12100, 52600], '350.6': [1214, 3950, 13200, 52600], '351.5': [1214, 3950, 13380, 52600], '352.4': [1214, 3950, 15500, 52700], '353.3': [1214, 3950, 16400, 52700], '354.2': [1214, 3950, 17000, 52700], '355.1': [1214, 3950, 17700, 52700], '356.0': [1214, 3950, 22700, 52700], '357.8': [1214, 3950, 23600, 52700]})


files = np.loadtxt('../scst/scst_list_firstscan.txt', dtype='str')

lowreg = []

for area in range(len(regions)):
    curregion = regions[area]
    print curregion

    t1 = Table.read('../Dunmaps/fwhm/scans/starcat_'+curregion+'mapweight_fwhm.txt', format='ascii')
    scan = Table(fits.open('../scst/'+files[area])[1].data)

    # Make list of all applicable regions
    lowreg.append(curregion)

    # Convert coordinates to gl, gb
    scangal = SkyCoord(scan['ra_acs']*u.deg, scan['dec_acs']*u.deg, frame='icrs').galactic
    med = np.median(scangal.l.degree)

    scancut = np.where((scangal.l.degree > med+0.3) | (scangal.l.degree < med-0.3) | (scan['gyro_rate_rss'] > 0.001))
    print scancut
    scan.remove_rows(scancut)

    scangal = SkyCoord(scan['ra_acs']*u.deg, scan['dec_acs']*u.deg, frame='icrs').galactic

    init = []
    end = []

    for reg in range(len(scan)):
        # Get first indices
        if (scan['NDCFEC'][reg] < 200000) & (len(init) == 0):
            init.append(reg)
            continue

        try:
            # If not at end of scan and HFEC region encountered, add that last index
            if (scan['NDCFEC'][reg] < 200000) & (scan['NDCFEC'][reg+1] > 200000) & (len(end) == 0):
                end.append(reg)
                continue

        except IndexError:
            # In the case that it reaches the end with no HFEC regions, add last index
            if (scan['NDCFEC'][reg] < 200000) & (len(end) == 0):
                end.append(reg)
                continue
            else:
                pass

        # If it enters a high FEC region, keep passing along
        if scan['NDCFEC'][reg] > 200000:
            continue

        # Now add regions in middle of scans that are right after HFEC regions
        if (scan['NDCFEC'][reg] < 200000) & (scan['NDCFEC'][reg-1] > 200000) & (len(init) != 0):
                init.append(reg)

        try:
            # If not at end of scan and HFEC region encountered, add that last index
            if (scan['NDCFEC'][reg] < 200000) & (scan['NDCFEC'][reg+1] > 200000) & (len(end) != 0):
                end.append(reg)

        except IndexError:
            # In the case that it reaches the end with no HFEC regions, add last index
            if (scan['NDCFEC'][reg] < 200000) & (len(end) != 0):
                end.append(reg)
            else:
                pass

    #print 'len(init) = ', len(init)
    #print 'len(end) = ', len(end)

    if len(init) != len(end):
        end.append(init[-1])

    table = Table()

    for ii in range(len(init)):
        indices = np.where((t1['gb'] > scangal.b.degree[end[ii]]) & (t1['gb'] < scangal.b.degree[init[ii]]))
        if len(indices[0]) == 0:
            indices = np.where((t1['gb'] < scangal.b.degree[end[ii]]) & (t1['gb'] > scangal.b.degree[init[ii]]))

        table = vstack([table, t1[indices]])

    ascii.write(table, '../lowfec/'+curregion+'_lowfec.txt', format='basic')

