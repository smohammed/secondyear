import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20
from astropy.io import fits, ascii

fourfields = 0
twofields = 1
onefield = 0

# Plot NUV comparisons
skyrange = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

mean0, median0, stdev0 = [], [], []
mean1, median1, stdev1 = [], [], []
mean2, median2, stdev2 = [], [], []
mean3, median3, stdev3 = [], [], []

alldata = fits.open('../Dunmaps/sex_galex_matches_total_mapweight.fits')[1].data

for region in range(0,360,5):
    if onefield == 1:
        a = alldata[np.where((alldata['gl_galex'] > region) & (alldata['gl_galex'] < region + 5))]
        m0, med0, std0 = [], [], []

        for magrange in np.arange(11.5, 22, 0.5):
            mag0 = np.where((a['nuv_galex'] > magrange) & (a['nuv_galex'] < magrange+1))
            m0.append(np.mean((a['nuv_sex']-a['nuv_galex'])[mag0]))
            std0.append(np.std((a['nuv_sex']-a['nuv_galex'])[mag0]))
            med0.append(np.median((a['nuv_sex']-a['nuv_galex'])[mag0]))

            mean0.append(np.mean((a['nuv_sex']-a['nuv_galex'])[mag0]))
            median0.append(np.median((a['nuv_sex']-a['nuv_galex'])[mag0]))
            stdev0.append(np.std((a['nuv_sex']-a['nuv_galex'])[mag0]))

        plt.scatter(a['nuv_galex'], a['nuv_sex']-a['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        plt.errorbar(np.arange(11.5, 22, 0.5), m0, yerr=std0, color='red', linewidth='2')
        plt.axhline(y=0, c='green')

        plt.xlabel('NUV$_{GAIS}$')
        plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        plt.title('gl ='+str(region)+'-'+str(region+5))
        plt.xlim((12, 23))
        plt.ylim((-3, 2))
        plt.xlim((12, 23))
        plt.ylim((-3, 2))

        #plt.show()
        plt.savefig('04-28-nuvcomp_sextractor_gais_'+str(region)+'.png')
        plt.clf()

    if twofields == 1:
        a = alldata[np.where((alldata['gl_galex'] > region) & (alldata['gl_galex'] < region + 5))]
        a0 = a[np.where((a['gb_galex'] > -10) & (a['gb_galex'] < 0))]
        a1 = a[np.where((a['gb_galex'] > 0) & (a['gb_galex'] < 10))]

        m0, med0, std0 = [], [], []
        m1, med1, std1 = [], [], []

        for magrange in np.arange(11.5, 22, 0.5):
            mag0 = np.where((a0['nuv_galex'] > magrange) & (a0['nuv_galex'] < magrange+1))
            m0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            std0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            med0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

            mag1 = np.where((a1['nuv_galex'] > magrange) & (a1['nuv_galex'] < magrange+1))
            m1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            std1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            med1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

            mean0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            median0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            stdev0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

            mean1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            median1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            stdev1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

        fig, (ax0, ax1) = plt.subplots(2, sharex=True, sharey=True)

        ax0.scatter(a0['nuv_galex'], a0['nuv_sex']-a0['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax0.errorbar(np.arange(11.5, 22, 0.5), m0, yerr=std0, color='red', linewidth='2')
        ax0.axhline(y=0, c='green')
        ax1.scatter(a1['nuv_galex'], a1['nuv_sex']-a1['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax1.errorbar(np.arange(11.5, 22, 0.5), m1, yerr=std1, color='red', linewidth='2')
        ax1.axhline(y=0, c='green')

        ax1.set_xlabel('NUV$_{GAIS}$')
        ax0.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('gl ='+str(region)+'-'+str(region+5))
        ax0.set_xlim((12, 23))
        ax0.set_ylim((-3, 2))
        ax1.set_xlim((12, 23))
        ax1.set_ylim((-3, 2))

        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        #plt.show()
        plt.savefig('04-28-nuvcomp_sextractor_gl2split_gais_'+str(region)+'.png')
        plt.clf()

        '''
        w = np.array(np.arange(11.5, 22, 0.5).tolist()*len(range(0,360,5)))

        mean0 = np.array(mean0)
        median0 = np.array(median0)
        stdev0 = np.array(stdev0)
        mean0[np.isnan(mean0)] = 0
        median0[np.isnan(median0)] = 0
        stdev0[np.isnan(stdev0)] = 0
        mean1 = np.array(mean1)
        median1 = np.array(median1)
        stdev1 = np.array(stdev1)
        mean1[np.isnan(mean1)] = 0
        median1[np.isnan(median1)] = 0
        stdev1[np.isnan(stdev1)] = 0

        fig, (ax0, ax1) = plt.subplots(2, sharex=True, sharey=True)
        for i in range(0, 9):
            ax0.errorbar(w[i*21:i*21+21], mean0[i*21:i*21+21], yerr=stdev0[i*21:i*21+21])
            ax1.errorbar(w[i*21:i*21+21], mean1[i*21:i*21+21], yerr=stdev1[i*21:i*21+21])
        ax0.axhline(y=0, c='green')
        ax1.axhline(y=0, c='green')

        ax1.set_xlabel('NUV$_{GAIS}$')
        ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('All available regions')
        ax0.set_xlim((11, 22))
        ax0.set_ylim((-2, 1))
        ax1.set_xlim((11, 22))
        ax1.set_ylim((-2, 1))
        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        plt.savefig('../images/03-10-nuvcomp_sextractor_gl2split_gais_mean.png')
        '''
    if fourfields == 1:
        a = alldata[np.where((alldata['gl_galex'] > region) & (alldata['gl_galex'] < region + 5))]
        a0 = a[np.where((a['gb_galex'] > -10) & (a['gb_galex'] < -5))]
        a1 = a[np.where((a['gb_galex'] > -5) & (a['gb_galex'] < 0))]
        a2 = a[np.where((a['gb_galex'] > 0) & (a['gb_galex'] < 5))]
        a3 = a[np.where((a['gb_galex'] > 5) & (a['gb_galex'] < 10))]
        m0, med0, std0 = [], [], []
        m1, med1, std1 = [], [], []
        m2, med2, std2 = [], [], []
        m3, med3, std3 = [], [], []

        for magrange in np.arange(11.5, 22, 0.5):
            mag0 = np.where((a0['nuv_galex'] > magrange) & (a0['nuv_galex'] < magrange+0.5))
            m0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            std0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            med0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            mag1 = np.where((a1['nuv_galex'] > magrange) & (a1['nuv_galex'] < magrange+0.5))
            m1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            std1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            med1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            mag2 = np.where((a2['nuv_galex'] > magrange) & (a2['nuv_galex'] < magrange+0.5))
            m2.append(np.mean((a2['nuv_sex']-a2['nuv_galex'])[mag2]))
            std2.append(np.std((a2['nuv_sex']-a2['nuv_galex'])[mag2]))
            med2.append(np.median((a2['nuv_sex']-a2['nuv_galex'])[mag2]))
            mag3 = np.where((a3['nuv_galex'] > magrange) & (a3['nuv_galex'] < magrange+0.5))
            m3.append(np.mean((a3['nuv_sex']-a3['nuv_galex'])[mag3]))
            std3.append(np.std((a3['nuv_sex']-a3['nuv_galex'])[mag3]))
            med3.append(np.median((a3['nuv_sex']-a3['nuv_galex'])[mag3]))

            mean0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            median0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            stdev0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

            mean1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            median1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            stdev1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

            mean2.append(np.mean((a2['nuv_sex']-a2['nuv_galex'])[mag2]))
            median2.append(np.median((a2['nuv_sex']-a2['nuv_galex'])[mag2]))
            stdev2.append(np.std((a2['nuv_sex']-a2['nuv_galex'])[mag2]))

            mean3.append(np.mean((a3['nuv_sex']-a3['nuv_galex'])[mag3]))
            median3.append(np.median((a3['nuv_sex']-a3['nuv_galex'])[mag3]))
            stdev3.append(np.std((a3['nuv_sex']-a3['nuv_galex'])[mag3]))

        fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)
        ax0.scatter(a0['nuv_galex'], a0['nuv_sex']-a0['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax0.errorbar(np.arange(11.5, 22, 0.5), m0, yerr=std0, color='red', linewidth='2')
        ax0.axhline(y=0, c='green')
        ax1.scatter(a1['nuv_galex'], a1['nuv_sex']-a1['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax1.errorbar(np.arange(11.5, 22, 0.5), m1, yerr=std1, color='red', linewidth='2')
        ax1.axhline(y=0, c='green')
        ax2.scatter(a2['nuv_galex'], a2['nuv_sex']-a2['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax2.errorbar(np.arange(11.5, 22, 0.5), m2, yerr=std2, color='red', linewidth='2')
        ax2.axhline(y=0, c='green')
        ax3.scatter(a3['nuv_galex'], a3['nuv_sex']-a3['nuv_galex'], alpha=0.1, edgecolor='none', facecolor='black')
        ax3.errorbar(np.arange(11.5, 22, 0.5), m3, yerr=std3, color='red', linewidth='2')
        ax3.axhline(y=0, c='green')

        ax3.set_xlabel('NUV$_{GAIS}$')
        ax2.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('gl ='+region)
        ax0.set_xlim((11, 22))
        ax0.set_ylim((-2, 1))
        ax1.set_xlim((11, 22))
        ax1.set_ylim((-2, 1))
        ax2.set_xlim((11, 22))
        ax2.set_ylim((-2, 1))
        ax3.set_xlim((11, 22))
        ax3.set_ylim((-2, 1))
        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        plt.savefig('../images/03-10-nuvcomp_sextractor_gl4split_gais_'+region+'.png')
        plt.clf()

        w = np.array(np.arange(11.5, 22, 0.5).tolist()*len(range(0,360,5)))
        mean0 = np.array(mean0)
        median0 = np.array(median0)
        stdev0 = np.array(stdev0)
        mean0[np.isnan(mean0)] = 0
        median0[np.isnan(median0)] = 0
        stdev0[np.isnan(stdev0)] = 0
        mean1 = np.array(mean1)
        median1 = np.array(median1)
        stdev1 = np.array(stdev1)
        mean1[np.isnan(mean1)] = 0
        median1[np.isnan(median1)] = 0
        stdev1[np.isnan(stdev1)] = 0
        mean2 = np.array(mean2)
        median2 = np.array(median2)
        stdev2 = np.array(stdev2)
        mean2[np.isnan(mean2)] = 0
        median2[np.isnan(median2)] = 0
        stdev2[np.isnan(stdev2)] = 0
        mean3 = np.array(mean3)
        median3 = np.array(median3)
        stdev3 = np.array(stdev3)
        mean3[np.isnan(mean3)] = 0
        median3[np.isnan(median3)] = 0
        stdev3[np.isnan(stdev3)] = 0

        fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)
        for i in range(0, 9):
            ax0.errorbar(w[i*21:i*21+21], mean0[i*21:i*21+21], yerr=stdev0[i*21:i*21+21])
            ax1.errorbar(w[i*21:i*21+21], mean1[i*21:i*21+21], yerr=stdev1[i*21:i*21+21])
            ax2.errorbar(w[i*21:i*21+21], mean2[i*21:i*21+21], yerr=stdev2[i*21:i*21+21])
            ax3.errorbar(w[i*21:i*21+21], mean3[i*21:i*21+21], yerr=stdev3[i*21:i*21+21])
        ax0.axhline(y=0, c='green')
        ax1.axhline(y=0, c='green')
        ax2.axhline(y=0, c='green')
        ax3.axhline(y=0, c='green')

        ax3.set_xlabel('NUV$_{GAIS}$')
        ax2.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('All available regions')
        ax0.set_xlim((11, 22))
        ax0.set_ylim((-2, 1))
        ax1.set_xlim((11, 22))
        ax1.set_ylim((-2, 1))
        ax2.set_xlim((11, 22))
        ax2.set_ylim((-2, 1))
        ax3.set_xlim((11, 22))
        ax3.set_ylim((-2, 1))

        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        plt.savefig('../images/02-18-nuvcomp_sextractor_gl4split_gais_mean.png')