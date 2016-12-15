from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 17

cat = 0
gais = 1

ext = 1

############################################################
# APOKASK and Bovy with gaia CMD
############################################################
if cat == 1:
    sg = fits.open('../sex_gaia_dust_interp.fits')[1].data
    sggal = SkyCoord(sg['gl_sex']*u.deg, sg['gb_sex']*u.deg, frame='galactic')

if gais == 1:
    sg = fits.open('../gais_gaia_dust.fits')[1].data
    sggal = SkyCoord(sg['gl_gais']*u.deg, sg['gb_gais']*u.deg, frame='galactic')


apo = fits.open('../APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('../BovyRC_TGAS.fits')[1].data
bovgal = SkyCoord(bov['l']*u.deg, bov['b']*u.deg, frame='galactic')
sg1ind, apoind, angsep1, ang3d = search_around_sky(sggal, apogal, 3*u.arcsec)
sg2ind, bovind, angsep2, ang3d = search_around_sky(sggal, bovgal, 3*u.arcsec)
c1 = hstack([Table(sg)[sg1ind], Table(apo)[apoind]])
c2 = hstack([Table(sg)[sg2ind], Table(bov)[bovind]])
c1['angsep'] = angsep1
c2['angsep'] = angsep2

cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')
agerange = np.unique(cmd['logage'])
colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

if ext == 1:
    scatter_contour((sg['nuv_mag']-sg['ebv']*7.76)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['Mg']-sg['ebv']*3.303), threshold=2000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.3), contour_args=dict(cmap=cm.gray))

    plt.scatter((c2['nuv_mag']-c2['ebv']*7.76)-(c2['phot_g_mean_mag_1']-c2['ebv']*3.303), c2['Mg']-c2['ebv']*3.303, s=80, label='Bovy', c=c2['C_H']/c2['FE_H'],  marker='s')
    plt.scatter((c1['nuv_mag']-c1['ebv']*7.76)-(c1['phot_g_mean_mag_1']-c1['ebv']*3.303), c1['Mg']-c1['ebv']*3.303, s=80, label='APO', c=c1['C_FE'])

if ext == 0:
    scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['Mg'], threshold=2000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.3), contour_args=dict(cmap=cm.gray))

    plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag_1'], c2['Mg'], edgecolor='none', s=80, label='Bovy', c=c2['C_H']/c2['FE_H'], marker='s')
    plt.scatter(c1['nuv_mag']-c1['phot_g_mean_mag_1'], c1['Mg'], edgecolor='none', s=80, label='APO', c=c1['C_FE'])

'''
for age in range(len(agerange)):
    cmd2 = cmd[np.where(cmd['logage'] == agerange[age])]
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
'''

# Just for a general plot
'''
plt.scatter((c2['nuv_mag']-c2['ebv']*7.76)-(c2['phot_g_mean_mag_1']-c2['ebv']*3.303), c2['Mg']-c2['ebv']*3.303, edgecolor='none', s=80, label='RC', c='red')
plt.scatter((c1['nuv_mag']-c1['ebv']*7.76)-(c1['phot_g_mean_mag_1']-c1['ebv']*3.303), c1['Mg']-c1['ebv']*3.303, edgecolor='none', s=80, c='red')

plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag_1'], c2['Mg'], edgecolor='none', s=80, label='RC', c='red')
plt.scatter(c1['nuv_mag']-c1['phot_g_mean_mag_1'], c1['Mg'], edgecolor='none', s=80, c='red')
'''

plt.xlabel('(NUV - E$_{B-V}$ * 7.76) - (G - E$_{B-V}$ * 3.303)')
plt.ylabel('MG - E$_{B-V}$ * 3.303')
plt.xlabel('NUV - G')
plt.ylabel('MG')
plt.colorbar().set_label('C/Fe')

if ext == 1:
    if cat == 1:
        plt.title('S+G, GSF15 ext with RC stars')
    if gais == 1:
        plt.title('GAIS+G, GSF15 ext with RC stars')

if ext == 0:
    if cat == 1:
        plt.title('S+G with RC stars')
    if gais == 1:
        plt.title('GAIS+G with RC stars')

plt.legend(scatterpoints=1, loc=3)
plt.xlim((-1, 11))
plt.ylim((8, -6))
plt.show()
