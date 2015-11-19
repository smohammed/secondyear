# Plot WD map

from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table, hstack, vstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

wds = Table.read('GAIS_vphas_2mass_wd.txt', format='ascii')

glval, gbval, numwd = [], [], []

for i in range(0, 360):
    for j in range(-10, 10):
        if len(np.where((wds['gl_galex'] > i) & (wds['gl_galex'] < i+1) & (wds['gb_galex'] > j) & (wds['gb_galex'] < j+1))[0]) > 0.:
            glval.append(i)
            gbval.append(j)
            numwd.append(len(np.where((wds['gl_galex'] > i) & (wds['gl_galex'] < i+1) & (wds['gb_galex'] > j) & (wds['gb_galex'] < j+1))[0]))

glval = np.array(glval)
gbval = np.array(gbval)
numwd = np.array(numwd)


wds = Table.read('GAIS_vphas_2mass_wd.txt', format='ascii')
gl = wds['gl_galex'].tolist()
gb = wds['gb_galex'].tolist()
dresult = query(gl,gb,coordsys='gal')

dist = 10**(1.+np.array(dresult['distmod'])/5.)/1000.


for i in range(len(dresult['best'][0])):    # 31 different distances
    avnum = []
    for j in range(len(wds)):               # 3525 white dwarfs
        avnum.append(dresult['best'][j][i])
    avnum = np.array(avnum) * 9.
    
    fig,axes = plt.subplots(2,2,sharey=True)
    axes[0,0].scatter(wds['gl_galex'],avnum,c='red',label='All WDs')
    axes[0,1].scatter(wds['gl_galex'],avnum,c='red')
    axes[1,1].scatter(wds['gl_galex'],avnum,c='red')
    axes[1,0].scatter(wds['gl_galex'],avnum,c='red')
    
    axes[0,0].set_xlim(-2,40)
    axes[0,1].set_xlim(205,252)
    axes[1,0].set_xlim(252,300)
    axes[1,1].set_xlim(300,362)
    '''axes[0,0].set_ylim(-6,6)
    axes[0,1].set_ylim(-6,6)
    axes[1,0].set_ylim(-6,6)
    axes[1,1].set_ylim(-6,6)'''
    axes[1,0].set_xlabel('gl')
    axes[1,1].set_xlabel('gl')
    axes[1,0].set_ylabel('A$_V$')
    axes[0,0].set_ylabel('A$_V$')
    axes[0,0].legend(scatterpoints=1)
    fig.subplots_adjust(wspace=0)
    axes[0,0].set_title('Distance = '+str(dist[i])+' kpc')
    plt.savefig('wdplots/wds_d'+str(dist[i])[:4])+'kpc.png')
    #plt.show()
    
