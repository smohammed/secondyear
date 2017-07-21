from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 14, 12
matplotlib.rcParams['font.size'] = 18
from colormaps import *

cat = 0
gais = 1

distcut = 0

cbarax = 'FE_H'
#cbarax = 'lnM'
#cbarax = 'lnAge'

if cbarax == 'FE_H':
    vmin = -0.5
    vmax = .35

if cbarax == 'lnM':
    vmin = -0.5
    vmax = 1.25

if cbarax == 'lnAge':
    vmin = -1
    vmax = 3

############################################################
# APOKASK and Bovy with gaia CMD
############################################################
if cat == 1:
    sg = fits.open('../sex_gaia_dust_interp.fits')[1].data
    sggal = SkyCoord(sg['gl_sex']*u.deg, sg['gb_sex']*u.deg, frame='galactic')

if gais == 1:
    sg = fits.open('../gais_tgas_apass_dust.fits')[1].data
    sg = sg[~np.isnan(sg['ebv'])]
    sggal = SkyCoord(sg['gl_gais']*u.deg, sg['gb_gais']*u.deg, frame='galactic')

'''
apo = fits.open('../APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('../BovyRC_TGAS.fits')[1].data
bovgal = SkyCoord(bov['l']*u.deg, bov['b']*u.deg, frame='galactic')
sg1ind, apoind, angsep1, ang3d = search_around_sky(sggal, apogal, 3*u.arcsec)
sg2ind, bovind, angsep2, ang3d = search_around_sky(sggal, bovgal, 3*u.arcsec)
c1 = hstack([Table(sg)[sg1ind], Table(apo)[apoind]])
rc = hstack([Table(sg)[sg2ind], Table(bov)[bovind]])
c1['angsep'] = angsep1
rc['angsep'] = angsep2
c1.remove_column('phot_g_mean_mag_2')
c1.rename_column('phot_g_mean_mag_1', 'phot_g_mean_mag')
rc.remove_column('phot_g_mean_mag_2')
rc.rename_column('phot_g_mean_mag_1', 'phot_g_mean_mag')
#rc = fits.open('../gais_tgas_dust_ness.fits')[1].data
'''
rc = Table.read('../rcall_match.txt', format='ascii')

#cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')
#agerange = np.unique(cmd['logage'])
#colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

if distcut == 1:
    xlow = np.sort(rc['dist'])[0]
    xhi = np.sort(rc['dist'])[-1]
    sg = sg[np.where((sg['dist'] > xlow) & (sg['dist'] < xhi))]

    #nuvgbins = 28
    nuvgbins = 36
    nuvgth = 1000
    bvbins = 80
    bvth = 1000


if distcut == 0:
    nuvgbins = 40
    nuvgth = 1000
    bvbins = 80
    bvth = 1000


fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)
# NUV - G with ext
#scatter_contour((sg['nuv_mag']-sg['ebv']*7.76)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['Mg']-sg['ebv']*3.303), threshold=nuvgth, log_counts=True, histogram2d_args=dict(bins=nuvgbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[1, 1])

axes[1, 1].scatter((sg['nuv_mag']-sg['ebv']*7.24)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['MG']-sg['ebv']*3.303), edgecolor='none', color='k', s=1, alpha=0.1)

cbar = axes[1, 1].scatter((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303), rc['MG']-rc['ebv']*3.303, s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

# NUV - G with no ext
#scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['Mg'], threshold=nuvgth, log_counts=True, histogram2d_args=dict(bins=nuvgbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[0, 1])

axes[0, 1].scatter(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['MG'], edgecolor='none', color='k', s=1, alpha=0.1)


axes[0, 1].scatter(rc['nuv_mag']-rc['phot_g_mean_mag'], rc['MG'], s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)


# B - V with ext
#scatter_contour((sg['B_AB']-sg['ebv']*3.626)-(sg['V_AB']-sg['ebv']*2.742), (sg['MV']-sg['ebv']*2.742), threshold=bvth, log_counts=True, histogram2d_args=dict(bins=bvbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[1, 0])

axes[1,0].scatter((sg['B_apass']-sg['ebv']*3.626)-(sg['V_apass']-sg['ebv']*2.742), (sg['V_apass']-sg['distmod']-sg['ebv']*2.742), edgecolor='none', color='k', s=1, alpha=0.1)

axes[1, 0].scatter((rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742), rc['V_apass']-rc['distmod']-rc['ebv']*2.742, s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

# B - V with no ext
#scatter_contour(sg['B_AB']-sg['V_AB'], sg['MV'], threshold=bvth, log_counts=True, histogram2d_args=dict(bins=bvbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[0, 0])

axes[0,0].scatter(sg['B_apass']-sg['V_apass'], sg['V_apass']-sg['distmod'], edgecolor='none', color='k', s=1, alpha=0.1)

axes[0, 0].scatter(rc['B_apass']-rc['V_apass'], rc['V_apass']-rc['distmod'], s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

'''
axes[1, 0].set_xlabel('(B - E$_{B-V}$ * 3.626) - (V - E$_{B-V}$ * 2.742)')
axes[1, 1].set_xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
axes[1, 0].set_ylabel('M$_{\lambda}$ - E$_{B-V}$ * A$_{\lambda}$')
axes[0, 0].set_ylabel('M$_{\lambda}$')
axes[0, 0].text(-0.4, 7.8, '$\lambda$ = V, E$_{B-V}$ = 0')
axes[0, 1].text(1.1, 7.8, '$\lambda$ = G, E$_{B-V}$ = 0')
axes[1, 0].text(-0.4, 7.8, '$\lambda$ = V, A$_{\lambda}$ = 2.742')
axes[1, 1].text(1.1, 7.8, '$\lambda$ = G, A$_{\lambda}$ = 3.303')
'''

axes[1, 0].set_xlabel('(B - V)$_0$')
axes[1, 1].set_xlabel('(NUV - G)$_0$')
axes[1, 0].set_ylabel('M$_{\lambda_0}$')
axes[0, 0].set_ylabel('M$_{\lambda}$')
axes[0, 0].text(-0.5, 7.8, '$\lambda$ = V, E$_{B-V}$ = 0')
axes[0, 1].text(1.1, 7.8, '$\lambda$ = G, E$_{B-V}$ = 0')
axes[1, 0].text(-0.5, 7.8, '$\lambda$ = V')
axes[1, 1].text(1.1, 7.8, '$\lambda$ = G')

'''
if distcut == 1:
    plt.suptitle('GAIS + TGAS, GSF extinction with Ness+16 RC stars, 170 < pc < 2100')

if distcut == 0:
    plt.suptitle('GAIS + TGAS, GSF extinction with Bovy+14 RC stars')
'''

fig.subplots_adjust(right=0.84)
cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])

if cbarax == 'FE_H':
    fig.colorbar(cbar, cax=cbar_ax, label='[Fe/H]')
if cbarax == 'lnM':
    fig.colorbar(cbar, cax=cbar_ax, label='lnM [Msol]')
if cbarax == 'lnAge':
    fig.colorbar(cbar, cax=cbar_ax, label='lnAge [Stellar age]')



#cm.set_label('Fe/H')
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
axes[0, 0].get_xaxis().set_ticklabels([])
axes[0, 0].get_yaxis().set_ticklabels(['', '6', '4', '2', '0', '-2'])
axes[0, 1].get_xaxis().set_ticklabels([])

axes[1, 1].set_xlim((1, 11.5))
axes[1, 1].set_ylim((7.99, -2))
axes[0, 1].set_xlim((1, 11.5))
axes[0, 1].set_ylim((7.99, -2))

axes[1, 0].set_xlim((-0.5, 2.0))
axes[1, 0].set_ylim((7.99, -2))
axes[0, 0].set_xlim((-0.5, 2.0))
axes[0, 0].set_ylim((7.99, -2))

plt.show()
