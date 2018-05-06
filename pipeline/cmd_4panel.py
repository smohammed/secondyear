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

if gais == 1:
    #sg = fits.open('../gais_tgas_apass_dust.fits')[1].data
    #sg = sg[~np.isnan(sg['ebv'])]
    asc = fits.open('../galex-asc-gaia-match_smaller.fits',memmap=True)[1].data
    pc = 1000./asc['parallax']
    negpar = np.where(pc > 0)
    pc = pc[negpar]
    nuv = asc['mag_nuv']
    g = asc['phot_g_mean_mag']
    bp = asc['phot_bp_mean_mag']
    rp = asc['phot_rp_mean_mag']
    ag = asc['a_g_val']
    pc = 1000./asc['parallax']
    negpar = np.where(pc > 0)
    pc = pc[negpar]
    nuv = nuv[negpar]
    g = g[negpar]
    bp = bp[negpar]
    rp = rp[negpar]
    ag = ag[negpar]
    distmod = 5. * np.log10(pc) - 5



rc = fits.open('../asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

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


fig, axes = plt.subplots(nrows=3, ncols=2)#, sharey=True)
# Bottom row will be plot with no RC stars
# panel e
axes[2, 0].scatter((bp-ag*3.626)-(rp-ag*2.742), (g-distmod-ag*2.742), edgecolor='none', color='k', s=1, alpha=0.1)
# panel f
axes[2, 1].scatter((nuv-ag*7.24)-(g-ag*3.303), (g-distmod-ag*3.303), edgecolor='none', color='k', s=1, alpha=0.1)

# B - V with no ext
# panel a
axes[0, 0].scatter(bp-rp, rp-distmod, edgecolor='none', color='k', s=1, alpha=0.1)
axes[0, 0].scatter(rc['phot_bp_mean_mag']-rc['phot_bp_mean_mag'], rc['phot_g_mean_mag']-rc['distmod'], s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

# B - V with ext
# panel c
axes[1, 0].scatter((bp-ag*3.626)-(rp-ag*2.742), (g-distmod-ag*2.742), edgecolor='none', color='k', s=1, alpha=0.1)
axes[1, 0].scatter((rc['phot_bp_mean_mag']-rc['ebv']*3.626)-(rc['phot_rp_mean_mag']-rc['ebv']*2.742), rc['phot_g_mean_mag']-rc['distmod']-rc['ebv']*2.742, s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

# NUV - G with no ext
# panel b
axes[0, 1].scatter(nuv-g, g-distmod, edgecolor='none', color='k', s=1, alpha=0.1)

axes[0, 1].scatter(rc['nuv_mag']-rc['phot_g_mean_mag'], rc['phot_g_mean_mag']-rc['distmod'], s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)

# NUV - G with ext
# panel d
axes[1, 1].scatter((nuv-ag*7.24)-(g-ag*3.303), (g-distmod-ag*3.303), edgecolor='none', color='k', s=1, alpha=0.1)

cbar = axes[1, 1].scatter((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303), rc['MG']-rc['ebv']*3.303, s=40, edgecolor='none', c=rc[cbarax], vmin=vmin, vmax=vmax, cmap=viridis)


axes[0, 0].set_ylabel('M$_{G}$')
axes[2, 0].set_xlabel('(G$_{BP}$ - G$_{RP}$)$_0$')
axes[1, 0].set_ylabel('M$_{G_0}$')
axes[2, 0].set_ylabel('M$_{G_0}$')
axes[2, 1].set_xlabel('(NUV - G)$_0$')

axes[0, 0].text(-0.5, 7.8, '$G$ = G, E$_{BP-RP}$ = 0')
axes[0, 1].text(1.1, 7.8, '$G$ = G, E$_{BP-RP}$ = 0')
axes[1, 0].text(-0.5, 7.8, '$G$ = G')
axes[1, 1].text(1.1, 7.8, '$G$ = G')
axes[2, 0].text(-0.5, 7.8, '$G$ = V')
axes[2, 1].text(1.1, 7.8, '$G$ = G')


# Label panels
axes[0, 0].text(1.7, 7.9, '(a)')
axes[0, 1].text(10.7, 7.85, '(b)')
axes[1, 0].text(1.7, 7.9, '(c)')
axes[1, 1].text(10.7, 7.85, '(d)')
axes[2, 0].text(1.7, 7.9, '(e)')
axes[2, 1].text(10.7, 7.9, '(f)')


fig.subplots_adjust(right=0.84)
cbar_ax = fig.add_axes([0.85, 0.15, 0.04, 0.7])

if cbarax == 'FE_H':
    fig.colorbar(cbar, cax=cbar_ax, label='[Fe/H]')
if cbarax == 'lnM':
    fig.colorbar(cbar, cax=cbar_ax, label='lnM [Msol]')
if cbarax == 'lnAge':
    fig.colorbar(cbar, cax=cbar_ax, label='lnAge [Stellar age]')

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
axes[0, 0].get_xaxis().set_ticklabels([])
axes[0, 1].get_xaxis().set_ticklabels([])
axes[1, 0].get_xaxis().set_ticklabels([])
axes[1, 1].get_xaxis().set_ticklabels([])

#axes[0, 0].get_yaxis().set_ticklabels([])
axes[0, 1].get_yaxis().set_ticklabels([])
#axes[1, 0].get_yaxis().set_ticklabels([])
axes[1, 1].get_yaxis().set_ticklabels([])
axes[2, 1].get_yaxis().set_ticklabels([])

axes[0, 0].get_yaxis().set_ticklabels(['-2', '0', '2', '4', '6', ''])
axes[1, 0].get_yaxis().set_ticklabels(['-2', '0', '2', '4', '6', ''])
axes[2, 0].get_yaxis().set_ticklabels(['-2', '0', '2', '4', '6', '8'])


axes[0, 0].set_xlim((-0.5, 1.99))
axes[0, 0].set_ylim((7.99, -2))

axes[1, 0].set_xlim((-0.5, 1.99))
axes[1, 0].set_ylim((7.99, -2))

axes[1, 1].set_xlim((1, 11.6))
axes[1, 1].set_ylim((7.99, -2))

axes[0, 1].set_xlim((1, 11.6))
axes[0, 1].set_ylim((7.99, -2))

axes[2, 0].set_xlim((-0.5, 1.99))
axes[2, 0].set_ylim((7.99, -2))

axes[2, 1].set_xlim((1, 11.6))
axes[2, 1].set_ylim((7.99, -2))

plt.show()
