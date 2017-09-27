from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 14, 12
matplotlib.rcParams['font.size'] = 23

cat = 1
gais = 0

ext = 0
nuvg = 1
bv = 0

if nuvg == 1:
    sgxval = 'nuv_mag'
    sgebv = 'ebv'
    sgyval = 'phot_g_mean_mag'
    avx = 7.24  # NUV
    avy = 3.303  # G
    yvalabsmag = 'MNUV'

if bv == 1:
    sgxval = 'B_AB'
    sgebv = 'ebv'
    sgyval = 'V_AB'
    avx = 3.626  # B
    avy = 2.742  # V
    yvalabsmag = 'MV'

############################################################
# APOKASK and Bovy with gaia CMD
############################################################
if cat == 1:
    sg = fits.open('../sex_gaia_dust_interp.fits')[1].data
    sggal = SkyCoord(sg['gl_sex']*u.deg, sg['gb_sex']*u.deg, frame='galactic')

if gais == 1:
    sg = fits.open('../gais_tgas_match_dust.fits')[1].data

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
c1.remove_column('phot_g_mean_mag_2')
c1.rename_column('phot_g_mean_mag_1', 'phot_g_mean_mag')
c2.remove_column('phot_g_mean_mag_2')
c2.rename_column('phot_g_mean_mag_1', 'phot_g_mean_mag')

#c2 = fits.open('../gais_tgas_dust_ness.fits')[1].data


cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')
agerange = np.unique(cmd['logage'])
colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

if nuvg == 1:
    bins = 40
    threshold = 1000

if bv == 1:
    bins = 80
    threshold = 1000

if ext == 1:
    #scatter_contour((sg[sgxval]-sg[sgebv]*avx)-(sg[sgyval]-sg[sgebv]*avy), (sg[sgyval]-sg['distmod']-sg[sgebv]*avy), threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins), plot_args=dict(color='k', markersize=1, alpha=1), contour_args=dict(cmap=cm.gray))

    #plt.scatter((c1[sgxval]-c1[sgebv]*avx)-(c1[sgyval]-c1[sgebv]*avy), c1[sgyval]-c1['distmod']-c1[sgebv]*avy, s=60, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
    #plt.scatter((c2[sgxval]-c2[sgebv]*avx)-(c2[sgyval]-c2[sgebv]*avy), c2[sgyval]-c2['distmod']-c2[sgebv]*avy, s=60, edgecolor='none', label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)

    scatter_contour((sg['nuv_mag']-sg['ebv']*7.24)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['nuv_mag']-sg['distmod']-sg['ebv']*7.24), threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins), plot_args=dict(color='k', markersize=1, alpha=1), contour_args=dict(cmap=cm.gray))

    plt.scatter((c1['nuv_mag']-c1['ebv']*7.24)-(c1['phot_g_mean_mag']-c1['ebv']*3.303), (c1['nuv_mag']-c1['distmod']-c1['ebv']*7.24), s=60, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
    plt.scatter((c2['nuv_mag']-c2['ebv']*7.24)-(c2['phot_g_mean_mag']-c2['ebv']*3.303), (c2['nuv_mag']-c2['distmod']-c2['ebv']*7.24), s=60, edgecolor='none', label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)

if ext == 0:
    #scatter_contour(sg[sgxval]-sg[sgyval], sg[sgyval]-sg['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins), plot_args=dict(color='k', markersize=1, alpha=1), contour_args=dict(cmap=cm.gray))

    #plt.scatter(c2[sgxval]-c2[sgyval], c2[sgyval]-c2['distmod'], edgecolor='none', s=60, label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)
    #plt.scatter(c1[sgxval]-c1[sgyval], c1[sgyval]-c1['distmod'], edgecolor='none', s=60, label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)

    scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['nuv_mag']-sg['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins), plot_args=dict(color='k', markersize=1, alpha=1), contour_args=dict(cmap=cm.gray))

    plt.scatter(c1['nuv_mag']-c1['phot_g_mean_mag'], c1['nuv_mag']-c1['distmod'], s=60, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
    plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['nuv_mag']-c2['distmod'], s=60, edgecolor='none', label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)

'''
for age in range(len(agerange)):
    cmd2 = cmd[np.where((cmd['logage'] == agerange[age]) & (cmd['Z'] == 0.0304))]
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
    plt.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
'''

if nuvg == 1:
    if ext == 1:
        plt.xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
        plt.ylabel('M$_{G}$ - E$_{B-V}$ * 3.303')

        if cat == 1:
            plt.title('S+G, GSF15 ext with RC stars')
        if gais == 1:
            plt.title('GAIS + TGAS, GSF extinction with RC stars')

    if ext == 0:
        plt.xlabel('NUV - G')
        plt.ylabel('M$_{NUV}$')

        if cat == 1:
            plt.title('S+G with RC stars')
        if gais == 1:
            plt.title('GAIS + TGAS with RC stars')

if bv == 1:
    if ext == 1:
        plt.xlabel('(B - E$_{B-V}$ * 3.626) - (V - E$_{B-V}$ * 2.742)')
        plt.ylabel('MV - E$_{B-V}$ * 2.742')

        if cat == 1:
            plt.title('UGPS + TGAS, GSF extinction with RC stars')
        if gais == 1:
            plt.title('GAIS + TGAS, GSF extinction with RC stars')

    if ext == 0:
        plt.xlabel('B - V ')
        plt.ylabel('MV')

        if cat == 1:
            plt.title('UGPS + TGAS with RC stars')
        if gais == 1:
            plt.title('GAIS + TGAS with RC stars')

cm = plt.colorbar()
cm.set_label('Fe/H')

#plt.legend(scatterpoints=1, loc=3)
if nuvg == 1:
    plt.xlim((1, 11.5))
    plt.ylim((8, -3))

if bv == 1:
    plt.xlim((-1, 2))
    plt.ylim((8, -3))

plt.show()
