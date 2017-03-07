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

cat = 0
gais = 1

ext = 1
nuvg = 0
bv = 1
distcut = 0

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

cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')
agerange = np.unique(cmd['logage'])
colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

if distcut == 1:
    xlow = np.sort(c2['dist'])[0]
    xhi = np.sort(c2['dist'])[-1]
    sg = sg[np.where((sg['dist'] > xlow) & (sg['dist'] < xhi))]

    nuvgbins = 28
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
scatter_contour((sg['nuv_mag']-sg['ebv']*7.76)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['MNUV']-sg['ebv']*7.76), threshold=nuvgth, log_counts=True, histogram2d_args=dict(bins=nuvgbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[1, 1])

cbar = axes[1, 1].scatter((c1['nuv_mag']-c1['ebv']*7.76)-(c1['phot_g_mean_mag']-c1['ebv']*3.303), c1['MNUV']-c1['ebv']*3.303, s=40, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
axes[1, 1].scatter((c2['nuv_mag']-c2['ebv']*7.76)-(c2['phot_g_mean_mag']-c2['ebv']*3.303), c2['MNUV']-c2['ebv']*3.303, s=40, edgecolor='none', label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)

# NUV - G with no ext
scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['MNUV'], threshold=nuvgth, log_counts=True, histogram2d_args=dict(bins=nuvgbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[0, 1])

axes[0, 1].scatter(c1['nuv_mag']-c1['phot_g_mean_mag'], c1['MNUV'], s=40, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
axes[0, 1].scatter(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['MNUV'], s=40, edgecolor='none', label='Bovy', c=c2['FE_H'], vmin=-0.5, vmax=0.35)


# B - V with ext
scatter_contour((sg['B_AB']-sg['ebv']*3.626)-(sg['V_AB']-sg['ebv']*2.742), (sg['MV']-sg['ebv']*2.742), threshold=bvth, log_counts=True, histogram2d_args=dict(bins=bvbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[1, 0])

axes[1, 0].scatter((c1['B_AB']-c1['ebv']*3.626)-(c1['V_AB']-c1['ebv']*2.742), c1['MV']-c1['ebv']*2.742, s=40, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
axes[1, 0].scatter((c2['B_AB']-c2['ebv']*3.626)-(c2['V_AB']-c2['ebv']*2.742), c2['MV']-c2['ebv']*2.742, s=40, edgecolor='none', label='APO', c=c2['FE_H'], vmin=-0.5, vmax=0.35)

# B - V with no ext
scatter_contour(sg['B_AB']-sg['V_AB'], sg['MV'], threshold=bvth, log_counts=True, histogram2d_args=dict(bins=bvbins), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=axes[0, 0])

axes[0, 0].scatter(c1['B_AB']-c1['V_AB'], c1['MV'], s=40, edgecolor='none', label='APO', c=c1['FE_H'], vmin=-0.5, vmax=0.35)
axes[0, 0].scatter(c2['B_AB']-c2['V_AB'], c2['MV'], s=40, edgecolor='none', label='APO', c=c2['FE_H'], vmin=-0.5, vmax=0.35)


axes[1, 0].set_xlabel('(B - E$_{B-V}$ * 3.626) - (V - E$_{B-V}$ * 2.742)')
axes[1, 1].set_xlabel('(NUV - E$_{B-V}$ * 7.76) - (NUV - E$_{B-V}$ * 7.76)')

axes[1, 0].set_ylabel('M$_{\lambda}$ - E$_{B-V}$ * A$_{\lambda}$')
axes[0, 0].set_ylabel('M$_{\lambda}$')

axes[0, 0].text(-0.9, 7.8, '$\lambda$ = V, E$_{B-V}$ = 0')
axes[0, 1].text(1.1, 7.8, '$\lambda$ = NUV, E$_{B-V}$ = 0')

axes[1, 0].text(-0.9, 7.8, '$\lambda$ = V, A$_{\lambda}$ = 2.742')
axes[1, 1].text(1.1, 7.8, '$\lambda$ = NUV, A$_{\lambda}$ = 3.303')

if distcut == 1:
    plt.suptitle('GAIS + TGAS, GSF extinction with RC stars, 380 < pc < 2000')

if distcut == 0:
    plt.suptitle('GAIS + TGAS, GSF extinction with RC stars')

fig.subplots_adjust(right=0.84)
cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])
fig.colorbar(cbar, cax=cbar_ax, label='Fe/H')


#cm.set_label('Fe/H')
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
axes[0, 0].get_xaxis().set_ticklabels([])
axes[0, 1].get_xaxis().set_ticklabels([])

axes[1, 1].set_xlim((1, 11.5))
axes[1, 1].set_ylim((7.99, -3))
axes[1, 1].set_xlim((1, 11.5))
axes[1, 1].set_ylim((7.99, -3))
axes[0, 1].set_xlim((1, 11.5))
axes[0, 1].set_ylim((7.99, -3))
axes[0, 1].set_xlim((1, 11.5))
axes[0, 1].set_ylim((7.99, -3))

axes[1, 0].set_xlim((-1, 1.99))
axes[1, 0].set_ylim((7.99, -3))
axes[1, 0].set_xlim((-1, 1.99))
axes[1, 0].set_ylim((7.99, -3))
axes[0, 0].set_xlim((-1, 1.99))
axes[0, 0].set_ylim((7.99, -3))
axes[0, 0].set_xlim((-1, 1.99))
axes[0, 0].set_ylim((7.99, -3))

plt.show()
