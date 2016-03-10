import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20
from astropy.io import fits, ascii

fourfields = 0
twofields = 1

# Plot NUV comparisons
skyrange = ['17.6-19.4', '20.3-25.7', '8.6-12.2', '205.7-210.2', '211.1-213.8', '214.7-217.4', '218.3-221.0', '223.7-226.4', '228.2-231.8']
mean0, median0, stdev0 = [], [], []
mean1, median1, stdev1 = [], [], []
mean2, median2, stdev2 = [], [], []
mean3, median3, stdev3 = [], [], []

alldata = fits.open('../sextractor_galex_matches_1-4.fits')[1].data

for region in range(0,360,5):
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
        ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('gl ='+str(region)+'-'+str(region+5))
        ax0.set_xlim((11, 22))
        ax0.set_ylim((-2, 1))
        ax1.set_xlim((11, 22))
        ax1.set_ylim((-2, 1))

        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        #plt.show()
        plt.savefig('../images/03-10-nuvcomp_sextractor_gl2split_gais_'+str(region)+'.png')
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
