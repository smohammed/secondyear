from astropy.io import fits
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 14, 12
matplotlib.rcParams['font.size'] = 18
from colormaps import *

cat = 0
gais = 1

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
    cat = fits.open('../sex_gaia_dust_interp.fits')[1].data

if gais == 1:
    cat = fits.open('../galex-asc-gaia-match_smaller.fits',memmap=True)[1].data
    pc = 1000./cat['parallax']
    negpar = np.where((pc > 0) & (cat['visibility_periods_used'] > 8) & (cat['parallax_error']/cat['parallax'] < 0.1) & (pc < 3500) & (cat['phot_bp_mean_mag'] > 0) & (cat['phot_rp_mean_mag'] > 0)) 
    nuv = cat['mag_nuv']
    g = cat['phot_g_mean_mag']
    bp = cat['phot_bp_mean_mag']
    rp = cat['phot_rp_mean_mag']
    ebv = cat['e_bv']

    pc = pc[negpar]
    nuv = nuv[negpar]
    g = g[negpar]
    bp = bp[negpar]
    rp = rp[negpar]
    ebv = ebv[negpar]

    pc = pc[~np.isnan(bp)]
    nuv = nuv[~np.isnan(bp)]
    g = g[~np.isnan(bp)]
    bp = bp[~np.isnan(bp)]
    rp = rp[~np.isnan(bp)]
    ebv = ebv[~np.isnan(bp)]
    pc = pc[~np.isnan(rp)]
    nuv = nuv[~np.isnan(rp)]
    g = g[~np.isnan(rp)]
    bp = bp[~np.isnan(rp)]
    rp = rp[~np.isnan(rp)]
    ebv = ebv[~np.isnan(rp)]

    distmod = 5. * np.log10(pc) - 5

rc = fits.open('../asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

#cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')
#agerange = np.unique(cmd['logage'])
#colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
# Panel a
scatter_contour(bp-rp, g-distmod, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-0.5,3], [-4,15]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 0])

# Panel b
scatter_contour((bp-ebv*2.85)-(rp-ebv*2.85), g-distmod-ebv*2.85, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-0.5,3], [-4,15]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 1])

# Panel c
scatter_contour(bp-rp, g-distmod, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-0.5,3], [-4,15]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 0])

axes[1, 0].scatter(rc['phot_bp_mean_mag']-rc['phot_rp_mean_mag'], rc['phot_g_mean_mag']-rc['distmod'], s=15, edgecolor='none', c=rc['Fe_H'], vmin=vmin, vmax=vmax, cmap=viridis, zorder=20)

# Panel d
scatter_contour((bp-ebv*2.85)-(rp-ebv*2.85), g-distmod-ebv*2.85, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-0.5,3], [-4,15]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 1])

cbar = axes[1, 1].scatter((rc['phot_bp_mean_mag']-rc['ebv']*2.85)-(rc['phot_rp_mean_mag']-rc['ebv']*2.85), rc['phot_g_mean_mag']-rc['distmod']-rc['ebv']*2.85, s=15, edgecolor='none', c=rc['Fe_H'], vmin=vmin, vmax=vmax, cmap=viridis, zorder=20)


axes[0, 0].set_xlim((-0.5, 3))
axes[0, 0].set_ylim((14, -3))

axes[0, 0].set_ylabel('M$_{G}$')
axes[1, 0].set_ylabel('M$_{G}$')

axes[1, 0].set_xlabel('(G$_{BP}$ - G$_{RP}$)')
axes[1, 1].set_xlabel('(G$_{BP}$ - G$_{RP}$)$_0$')

axes[0, 0].text(1.5, 13, 'E(B-V) = 0, (a)')
axes[0, 1].text(2.7, 13, '(b)')
axes[1, 0].text(1.5, 13, 'E(B-V) = 0, (c)')
axes[1, 1].text(2.7, 13, '(d)')


fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=0.84)
cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])
fig.colorbar(cbar, cax=cbar_ax, label='[Fe/H]')
plt.show()




# NUV - G
fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
# Panel a
scatter_contour(nuv-g, g-distmod, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 0])

# Panel b
scatter_contour((nuv-ebv*7.24)-(g-ebv*2.85), g-distmod-ebv*2.85, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 1])

# Panel c
scatter_contour(nuv-g, g-distmod, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 0])

axes[1, 0].scatter(rc['nuv_mag']-rc['phot_g_mean_mag'], rc['phot_g_mean_mag']-rc['distmod'], s=15, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis, zorder=10)

# Panel d
scatter_contour((nuv-ebv*7.24)-(g-ebv*2.85), g-distmod-ebv*2.85, threshold=1000, log_counts=True, histogram2d_args=dict(bins=100, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.01), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1,1])

cbar = axes[1, 1].scatter((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85), rc['MG']-rc['ebv']*2.85, s=15, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis, zorder=20)


axes[0, 0].set_xlim((-2, 11.6))
axes[0, 0].set_ylim((14, -3))

axes[1, 0].set_xlabel('(NUV - G)')
axes[1, 1].set_xlabel('(NUV - G)$_0$')
axes[0, 0].set_ylabel('M$_{G}$')
axes[1, 0].set_ylabel('M$_{G}$')


axes[0, 0].text(5, 13, 'E(B-V) = 0, (a)')
axes[0, 1].text(10, 13, '(b)')
axes[1, 0].text(5, 13, 'E(B-V) = 0, (c)')
axes[1, 1].text(10, 13, '(d)')


fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=0.84)
cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])
fig.colorbar(cbar, cax=cbar_ax, label='[Fe/H]')
plt.show()




