import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# Plot NUV comparisons
skyrange = ['17.6-19.4', '20.3-25.7', '8.6-12.2', '205.7-210.2', '211.1-213.8', '214.7-217.4', '218.3-221.0', '223.7-226.4', '228.2-231.8']
fourfields = 1
twofields = 0

mean = Table(np.zeros(21))
median = Table(np.zeros(21))
std = Table(np.zeros(21))

mean.remove_rows(0)
median.remove_rows(0)
std.remove_rows(0)

for region in skyrange:
    if fourfields == 1:
        a = Table.read('../sex_galex_matches_'+region.replace('.', '')+'.txt', format='ascii')
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

        mean = vstack([mean,Table(np.array(m0))])
        median = vstack([median,Table(np.array(m0))])
        std = vstack([std,Table(np.array(m0))])

        fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=True)
        ax0.scatter(a0['nuv_galex'], a0['nuv_sex']-a0['nuv_galex'], alpha=0.1,edgecolor='none',facecolor='black')
        ax0.errorbar(np.arange(11.5, 22,0.5), m0, yerr=std0, color='red', linewidth='2')
        ax0.axhline(y=0, c='green')
        ax1.scatter(a1['nuv_galex'], a1['nuv_sex']-a1['nuv_galex'], alpha=0.1,edgecolor='none',facecolor='black')
        ax1.errorbar(np.arange(11.5, 22,0.5), m1, yerr=std1, color='red', linewidth='2')
        ax1.axhline(y=0, c='green')
        ax2.scatter(a2['nuv_galex'], a2['nuv_sex']-a2['nuv_galex'], alpha=0.1,edgecolor='none',facecolor='black')
        ax2.errorbar(np.arange(11.5, 22,0.5), m2, yerr=std2, color='red', linewidth='2')
        ax2.axhline(y=0, c='green')
        ax3.scatter(a3['nuv_galex'], a3['nuv_sex']-a3['nuv_galex'], alpha=0.1,edgecolor='none',facecolor='black')
        ax3.errorbar(np.arange(11.5, 22,0.5), m3, yerr=std3, color='red', linewidth='2')
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

        plt.show()
        plt.savefig('../images/02-17-nuvcomp_sextractor_gl4split_gais_'+region+'.png')
        plt.clf()

    if twofields == 1:
        a = Table.read('../sex_galex_matches_'+region.replace('.', '')+'.txt', format='ascii')
        a0 = a[np.where((a['gb_galex'] > -10) & (a['gb_galex'] < 0))]
        a1 = a[np.where((a['gb_galex'] > 0) & (a['gb_galex'] < 10))]

        m0, med0, std0 = [], [], []
        m1, med1, std1 = [], [], []

        for magrange in np.arange(11.5, 22, 1):
            mag0 = np.where((a0['nuv_galex'] > magrange) & (a0['nuv_galex'] < magrange+1))
            m0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            std0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
            med0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

            mag1 = np.where((a1['nuv_galex'] > magrange) & (a1['nuv_galex'] < magrange+1))
            m1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            std1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
            med1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

        fig, (ax0, ax1) = plt.subplots(2, sharex=True, sharey=True)

        ax0.scatter(a0['nuv_galex'], a0['nuv_sex']-a0['nuv_galex'], alpha=0.1)
        ax0.errorbar(np.arange(11.5, 22), m0, yerr=std0, color='red', linewidth='2')
        ax0.axhline(y=0, c='green')
        ax1.scatter(a1['nuv_galex'], a1['nuv_sex']-a1['nuv_galex'], alpha=0.1)
        ax1.errorbar(np.arange(11.5, 22), m1, yerr=std1, color='red', linewidth='2')
        ax1.axhline(y=0, c='green')

        ax1.set_xlabel('NUV$_{GAIS}$')
        ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('gl ='+region)
        ax0.set_xlim((11, 22))
        ax0.set_ylim((-2, 1))
        ax1.set_xlim((11, 22))
        ax1.set_ylim((-2, 1))

        fig.subplots_adjust(hspace=0)
        plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
        #plt.show()
        plt.savefig('../images/02-17-nuvcomp_sextractor_gl2split_gais_'+region+'.png')
        plt.clf()
