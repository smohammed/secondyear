from astropy.io import fits
from astropy.table import Table
import matplotlib
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from astroML.plotting import scatter_contour
matplotlib.rcParams['figure.figsize'] = 14, 12
matplotlib.rcParams['font.size'] = 18
from colormaps import *

# Load original and matched catalogs
cat = fits.open('../starcat_allscans_10-12-18_cuts.fits')[1].data # catalog
ps = fits.open('../plane_ps1_g10-20_11_27_18.fits')[1].data # Pan starrs
cg = fits.open('../plane_gaiadr2_dust_11_14_18.fits')[1].data  # gaia
gg = fits.open('../plane_gais_11_27_18.fits')[1].data # GAIS 

############################################################
# Now make plots
############################################################
# NUV vs FWHM
scatter_contour(cat['nuv'], cat['FWHM_IMAGE'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10))
plt.xlim(13, 20)
plt.xlabel('NUV')
plt.ylabel('FWHM [pixels]')
plt.show()




############################################################
# Gaia CMD
############################################################
gag = fits.open('gais_gaiadr2_negparcuts_dust_09_27_18.fits')[1].data 

negpar = np.where((cg['dist'] > 0) & (cg['visibility_periods_used'] > 8) & (cg['phot_bp_mean_mag'] > 0) & (cg['phot_rp_mean_mag'] > 0) & (cg['expsum'] > 5) & (cg['ebv'] > 0) & (cg['parallax_error']/cg['parallax'] < 0.1))
#cut = np.where((gag['parallax'] > 0) & (gag['phot_g_mean_mag'] > 0) & (gag['phot_bp_mean_mag'] > 0) & (gag['phot_rp_mean_mag'] > 0) & (gag['visibility_periods_used'] > 8) & (gag['mag_nuv'] > 0) & (gag['glat'] > -10) & (gag['glat'] < 10))

cg = cg[negpar]
gag = gag[np.where((gag['ebv'] > 0) & (gag['parallax_error']/gag['parallax'] < 0.1))]

threshold = 1000
bins = 40

fig, axes = plt.subplots(nrows=2, ncols=2)

scatter_contour(gag['mag_nuv']-gag['phot_g_mean_mag'], gag['phot_g_mean_mag']-gag['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 0])

scatter_contour(cg['nuv']-cg['phot_g_mean_mag'], cg['phot_g_mean_mag']-cg['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 1])

scatter_contour((gag['mag_nuv']-gag['ebv']*7.24)-(gag['phot_g_mean_mag']-gag['ebv']*2.85), gag['phot_g_mean_mag']-gag['distmod']-gag['ebv']*2.85, threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 0])


scatter_contour((cg['nuv']-cg['ebv']*7.24)-(cg['phot_g_mean_mag']-cg['ebv']*2.85), cg['phot_g_mean_mag']-cg['distmod']-cg['ebv']*2.85, threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 1])


axes[1, 0].set_xlabel('NUV$_{GAIS}$ - G')
axes[1, 1].set_xlabel('NUV$_{Plane}$ - G')

axes[0, 0].set_ylabel('M$_G$')
axes[1, 0].set_ylabel('M$_{G_0}$')

axes[0, 0].text(8.4, 13.9, 'E(B-V) = 0')
axes[0, 1].text(8.4, 13.9, 'E(B-V) = 0')

axes[0, 0].set_xlim((-2, 11.6))
axes[0, 0].set_ylim((14, -3))
axes[1, 0].set_xlim((-2, 11.6))
axes[1, 0].set_ylim((14, -3))

axes[0, 1].set_xlim((-2, 11.6))
axes[0, 1].set_ylim((14, -3))
axes[1, 1].set_xlim((-2, 11.6))
axes[1, 1].set_ylim((14, -3))

axes[0, 0].set_xticks([])
axes[0, 1].set_yticks([])
axes[1, 1].set_yticks([])
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.show()


fig, axes = plt.subplots(nrows=2, ncols=2)
axes[0, 0].hist2d(gag['mag_nuv']-gag['phot_g_mean_mag'], gag['phot_g_mean_mag']-gag['distmod'], bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[0, 1].hist2d(cg['nuv']-cg['phot_g_mean_mag'], cg['phot_g_mean_mag']-cg['distmod'], bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[1, 0].hist2d((gag['mag_nuv']-gag['ebv']*7.24)-(gag['phot_g_mean_mag']-gag['ebv']*2.85), gag['phot_g_mean_mag']-gag['distmod']-gag['ebv']*2.85, bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[1, 1].hist2d((cg['nuv']-cg['ebv']*7.24)-(cg['phot_g_mean_mag']-cg['ebv']*2.85), cg['phot_g_mean_mag']-cg['distmod']-cg['ebv']*2.85, bins=(1000, 1000), norm=matplotlib.colors.LogNorm())


axes[1, 0].set_xlabel('NUV$_{GAIS}$ - G')
axes[1, 1].set_xlabel('NUV$_{Plane}$ - G')

axes[0, 0].set_ylabel('M$_G$')
axes[1, 0].set_ylabel('M$_{G_0}$')

axes[0, 0].text(8.4, 13.9, 'E(B-V) = 0')
axes[0, 1].text(8.4, 13.9, 'E(B-V) = 0')

axes[0, 0].set_xlim((-2, 11.6))
axes[0, 0].set_ylim((14, -5))
axes[1, 0].set_xlim((-2, 11.6))
axes[1, 0].set_ylim((14, -5))

axes[0, 1].set_xlim((-2, 11.6))
axes[0, 1].set_ylim((14, -5))
axes[1, 1].set_xlim((-2, 11.6))
axes[1, 1].set_ylim((14, -5))

axes[0, 0].set_xticks([])
axes[0, 1].set_yticks([])
axes[1, 1].set_yticks([])
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.show()


############################################################
# 4 panel grvsnuvg for paper
############################################################
from matplotlib.colors import LogNorm
ebvcut = np.where(psg['ebv']> 0)
nuvgmag = ((psg['nuv']-psg['ebv']*7.24)-(psg['gMeanPSFMag']-psg['ebv']*3.172))[ebvcut]
gr = ((psg['gMeanPSFMag']-psg['ebv']*3.172)-(psg['rMeanPSFMag']-psg['ebv']*2.271))[ebvcut]

pickles2 = pickles[np.arange(0, len(pickles), 5)]



fig, axes = plt.subplots(2, 2, sharex=True, sharey=True,figsize=(16, 12), facecolor='w')

cmap = axes[0, 0].hist2d(psg['nuv']-psg['gMeanPSFMag'], psg['gMeanPSFMag']-psg['rMeanPSFMag'], bins=5000,cmap='viridis', norm=LogNorm())
axes[0, 0].scatter(pickles2['nuv']-pickles2['g'], pickles2['g']-pickles2['r'], color='red', label='SED model',  s=30, zorder=20)

axes[0, 1].hist2d(psg['nuv']-psg['gMeanPSFMag'], psg['gMeanPSFMag']-psg['rMeanPSFMag'], bins=5000,cmap='viridis', norm=LogNorm())
axes[0, 1].scatter(wds['nuv']-wds['gMeanPSFMag'], wds['gMeanPSFMag']-wds['rMeanPSFMag'],s=1, label='WDs')

axes[1, 0].hist2d(psg['nuv']-psg['gMeanPSFMag'], psg['gMeanPSFMag']-psg['rMeanPSFMag'], bins=5000,cmap='viridis', norm=LogNorm())
axes[1, 0].scatter(ob['nuv']-ob['gMeanPSFMag'], ob['gMeanPSFMag']-ob['rMeanPSFMag'], s=1, label='Selected blue objects')

axes[1, 1].hist2d(nuvgmag, gr, bins=5000,cmap='viridis', norm=LogNorm())
axes[1, 1].scatter((ob['nuv']-ob['ebv']*7.24)-(ob['gMeanPSFMag']-ob['ebv']*3.172), (ob['gMeanPSFMag']-ob['ebv']*3.172)-(ob['rMeanPSFMag']-ob['ebv']*2.271), s=1)


axes[0, 0].arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
axes[0, 1].arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
axes[1, 0].arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
axes[1, 1].arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
axes[0, 0].annotate('O', xy=(0, -0.9), size=20)
axes[0, 0].annotate('B', xy=(1, -0.81), size=20)
axes[0, 0].annotate('A', xy=(2.5, -0.5), size=20)
axes[0, 0].annotate('F', xy=(3.7, -0.5), size=20)
axes[0, 0].annotate('G', xy=(5, -0.3), size=20)
axes[0, 0].annotate('K', xy=(6.2, -0.1), size=20)

axes[0, 0].legend(scatterpoints=1)
axes[0, 1].legend(scatterpoints=1)
axes[1, 0].legend(scatterpoints=1)

fig.subplots_adjust(wspace=0, hspace=0)
axes[0, 0].set_xlim((-1.2, 8))
axes[0, 0].set_ylim((-1.5, 2))
axes[1, 0].set_xlim((-1.2, 8))
axes[1, 0].set_ylim((-1.5, 2))
axes[1, 1].set_xlim((-1.2, 8))
axes[1, 1].set_ylim((-1.5, 2))
axes[1, 0].set_xlabel('NUV - g')
axes[1, 1].set_xlabel('NUV - g')
axes[0, 0].set_ylabel('g - r')
axes[1, 0].set_ylabel('g - r')


cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax)


plt.show()



############################################################
# Gaia NUVG vs MNUV
############################################################
fig, axes = plt.subplots(nrows=2, ncols=2)
axes[0, 0].hist2d(gag['mag_nuv']-gag['phot_g_mean_mag'], gag['mag_nuv']-gag['distmod'], bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[0, 1].hist2d(cg['nuv']-cg['phot_g_mean_mag'], cg['nuv']-cg['distmod'], bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[1, 0].hist2d((gag['mag_nuv']-gag['ebv']*7.24)-(gag['phot_g_mean_mag']-gag['ebv']*2.85), gag['mag_nuv']-gag['distmod']-gag['ebv']*7.24, bins=(1000, 1000), norm=matplotlib.colors.LogNorm())

axes[1, 1].hist2d((cg['nuv']-cg['ebv']*7.24)-(cg['phot_g_mean_mag']-cg['ebv']*2.85), cg['nuv']-cg['distmod']-cg['ebv']*7.24, bins=(1000, 1000), norm=matplotlib.colors.LogNorm())


axes[1, 0].set_xlabel('NUV$_{GAIS}$ - G')
axes[1, 1].set_xlabel('NUV$_{Plane}$ - G')

axes[0, 0].set_ylabel('M$_NUV$')
axes[1, 0].set_ylabel('M$_{NUV_0}$')

axes[0, 0].text(8.4, -3, 'E(B-V) = 0')
axes[0, 1].text(8.4, -3, 'E(B-V) = 0')

axes[0, 0].set_xlim((-2, 11.6))
axes[0, 0].set_ylim((18, -5))
axes[1, 0].set_xlim((-2, 11.6))
axes[1, 0].set_ylim((18, -5))

axes[0, 1].set_xlim((-2, 11.6))
axes[0, 1].set_ylim((18, -5))
axes[1, 1].set_xlim((-2, 11.6))
axes[1, 1].set_ylim((18, -5))

axes[0, 0].set_xticks([])
axes[0, 1].set_yticks([])
axes[1, 1].set_yticks([])
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.show()




############################################################
# CMD with Gaia and PS2
############################################################
cat = fits.open('../plane_gaiadr2_dust_06_12_19.fits')[1].data
cgp = fits.open('../plane_gaia_ps2_08_06_19.fits')[1].data


threshold = 1000
bins = 40

fig, axes = plt.subplots(nrows=2, ncols=2)

scatter_contour(cat['mag_nuv']-cat['phot_g_mean_mag'], cat['phot_g_mean_mag']-cat['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 0])

scatter_contour(cgp['nuv']-cgp['phot_g_mean_mag'], cgp['phot_g_mean_mag']-cgp['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 1])

scatter_contour((cat['mag_nuv']-cat['ebv']*7.24)-(cat['phot_g_mean_mag']-cat['ebv']*2.85), cat['phot_g_mean_mag']-cat['distmod']-cat['ebv']*2.85, threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 0])


scatter_contour((cgp['nuv']-cgp['ebv']*7.24)-(cgp['phot_g_mean_mag']-cgp['ebv']*2.85), cgp['phot_g_mean_mag']-cgp['distmod']-cgp['ebv']*2.85, threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 1])


axes[1, 0].set_xlabel('NUV$_{Gaia}$ - G')
axes[1, 1].set_xlabel('NUV$_{Gaia+PS2}$ - G')

axes[0, 0].set_ylabel('M$_G$')
axes[1, 0].set_ylabel('M$_{G_0}$')

axes[0, 0].text(8.4, 13.9, 'E(B-V) = 0')
axes[0, 1].text(8.4, 13.9, 'E(B-V) = 0')

axes[0, 0].set_xlim((-2, 11.6))
axes[0, 0].set_ylim((14, -3))
axes[1, 0].set_xlim((-2, 11.6))
axes[1, 0].set_ylim((14, -3))

axes[0, 1].set_xlim((-2, 11.6))
axes[0, 1].set_ylim((14, -3))
axes[1, 1].set_xlim((-2, 11.6))
axes[1, 1].set_ylim((14, -3))

axes[0, 0].set_xticks([])
axes[0, 1].set_yticks([])
axes[1, 1].set_yticks([])
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.show()


############################################################
# Angular separations for Gaia and PS1
############################################################
cg = fits.open('plane_gaiadr2_3ddust_ebvsfd_07_24_20.fits')[1].data
ps = fits.open('plane_ps2_08_01_19.fits')[1].data
ps = Table(ps)
dra = (ps['ALPHA_J2000']-ps['raMean'])*3600
ddec = (ps['DELTA_J2000']-ps['decMean'])*3600
angsep = np.sqrt(dra**2+ddec**2)
ps['angsep'] = angsep

cg = cg[np.where(cg['angsep']*3600. < 2.0)]
ps = ps[np.where(ps['angsep'] < 2.0)]

plt.hist(cg['angsep']*3600, bins=20, alpha=0.7, label='Gaia match', color='#00E6FF')
plt.hist(ps['angsep'], bins=20, alpha=0.7, label='PS2 match', color='#FF5C5C', histtype='step', linewidth=3, stacked=True)
plt.xlabel('Angular Separation')
plt.legend(scatterpoints=1)
plt.show()


############################################################
# NUV comparison
############################################################
gg = fits.open('plane_gais_06_12_19.fits')[1].data

# Subtract 0.25 mag from plane NUV because the offset was added to the _cut.fits file
scatter_contour(gg['nuv_mag'], gg['nuv_mag']-(gg['nuv']), threshold=10000, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))

plt.axhline(y=0, color='red')
plt.xlim(13, 21)
plt.ylim(-0.5, 1)
plt.xlabel('NUV (GAIS)')
plt.ylabel('$\Delta$NUV (GAIS - Plane)')
plt.show()


plt.hist2d(gg['nuv_mag'], gg['nuv_mag']-gg['nuv']-0.25, bins=(400,400))

plt.axhline(y=0.25, color='black', linestyle='--')
plt.axhline(y=0, color='red')
plt.xlim(13, 21)
plt.ylim(-0.5, 1)
plt.xlabel('NUV (GAIS)')
plt.ylabel('$\Delta$NUV (GAIS - Plane)')
plt.colorbar()
plt.show()

############################################################
# NUV histogram for all surveys
############################################################
cat = fits.open('starcat_allscans_06_12_19.fits')[1].data # catalog
ps = fits.open('plane_ps2_08_01_19.fits')[1].data # Pan starrs
cg = fits.open('plane_gaiadr2_dust_06_12_19.fits')[1].data  # gaia
gg = fits.open('plane_gais_08_26_19.fits')[1].data # GAIS 
allcat = fits.open('plane_gaiadr2_gais_ps2_09_04_19.fits')[1].data
none = Table.read('plane_innosurveys_08_16_19.txt',format='ascii')
#cg = cg[np.where(cg['angsep']*3600. < 2.0)]
#ps = ps[np.where((ps['angsep']*3600. < 2.0) & (ps['expsum'] > 3.) & (ps['ctsum'] > 1))]
#ps = ps[~np.isnan(cat['bkgdsum'])]

import matplotlib.ticker as ticker
fig, ax = plt.subplots()
ax.hist(cat['nuv'], bins=20, range=[12,21], label='UVGAPS', alpha=0.7, color='#00E6FF')
ax.hist(none['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='No surveys', color='yellow')
ax.hist(cg['nuv'], bins=20, range=[12, 21], histtype='step',linewidth=3, stacked=True, label='UVGAPS+Gaia')
ax.hist(gg['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='UVGAPS+GAIS')
ax.hist(ps['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=3, stacked=True, label='UVGAPS+PS2')
ax.hist(allcat['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='UVGAPS+All surveys')

ax.legend(loc=2)
ax.set_xlabel('NUV')
ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=len(cat)))
plt.show()

############################################################
# NUV offset histogram with GAIS
############################################################
gg = fits.open('plane_gais_08_26_19.fits')[1].data # GAIS 
nuv = gg['nuv']

g14 = (gg['nuv_mag']-nuv)[np.where((nuv > 16) & (nuv < 18))] - 0.25

plt.hist(g14, range=[-1, 1], bins=40)


for mag in range(12, 20,2):
    cut = np.where((nuv > mag) & (nuv < mag+2))
    dnuv = (gg['nuv_mag']-nuv)[cut] - 0.25
    plt.hist(dnuv, histtype='step', fill=False, stacked=True, label=str(mag)+'-'+str(mag+2), range=[-1, 1], bins=40)
plt.legend(scatterpoints=1, loc=2)
plt.xlabel('$\Delta$NUV (GAIS - UVGAPS)')
plt.xlim(-0.6, 0.6)
plt.show()


############################################################
# g vs g-r
############################################################
threshold = 10000
bins = 4

scatter_contour(ps['g_ps']-ps['r_ps'], ps['g_ps'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=(bins)), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.gca().invert_yaxis()
plt.xlim(-4, 4)
plt.ylim(20, 11)
plt.xlabel('g - r')
plt.ylabel('g')
plt.show()


############################################################
# g-r vs NUV-g
############################################################
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
fig, ax = plt.subplots()
pickles = Table.read('picklemags_laphare_final.txt', format='ascii')
scatter_contour(ps['nuv']-ps['gMeanPSFMag'], ps['gMeanPSFMag']-ps['rMeanPSFMag'], threshold=10000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.scatter(pickles['nuv']-pickles['g'], pickles['g']-pickles['r'], color='red', label='SED model',  s=30, zorder=10)
plt.xlim((-1.2, 7))
plt.ylim((-1, 1))

plt.annotate('O', xy=(0, -0.9), size=20)
plt.annotate('B', xy=(1, -0.81), size=20)
plt.annotate('A', xy=(2.5, -0.5), size=20)
plt.annotate('F', xy=(3.7, -0.5), size=20)
plt.annotate('G', xy=(5, -0.3), size=20)
plt.annotate('K', xy=(6.2, -0.1), size=20)
plt.xlabel('NUV - g')
plt.ylabel('g - r')
plt.arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
plt.legend(scatterpoints=1, loc=4)
plt.show()




# MS, WDs and binaries respectively
p1 = Polygon([[1.7, -0.5], [1.9, 0.3], [2.7, 0.7], [4, 1.05], [5.2, 1.05], [5.8, 1], [6.5, 0.95], [6.6, 0.3], [6, -0.1], [4.6, -0.4], [3.2, -0.5], [1.7, -0.5]])
p2 = Polygon([[-1, -0.6], [0.9, -0.25], [0.73, 0], [-1, -0.3], [-1, -0.6]])
p3 = Polygon([[1.7, -0.5], [1.9, 0.3], [2.7, 0.7], [4, 1.05], [5.2, 1.05], [5.8, 1], [6.5, 0.95], [6.5, 2], [-2.1, 2], [-2.1, 0.5], [1.9, 0.3]])

patches = [p1, p2, p3]
p = PatchCollection(patches, alpha=0.3, cmap=matplotlib.cm.jet)
p.set_array(np.array([100, 50, 80]))
#p.set_array(np.array(100.*np.random.rand(3)))
ax.add_collection(p)
plt.show()



# Or do this hist2d plot
fig, axes = plt.subplots(nrows=2, ncols=2)

axes[0, 0].hist2d(ps['nuv']-ps['gMeanPSFMag'], ps['gMeanPSFMag']-ps['rMeanPSFMag'], bins=(300,300),vmin=0, vmax=300)
plt.xlim(-2, 10)
plt.ylim(-1, 1.5)
plt.colorbar()
plt.show()


############################################################
# Area
############################################################
g = fits.open('gr67tile_table.fits')[1].data
ra = np.array(g['ra_cent'], dtype=np.float32)
dec = np.raray(g['dec_cent'], dtype=np.float32)
gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
gl = gal.galactic.l.degree
gb = gal.galactic.b.degree
q = np.where((gb > -10) & (gb < 10))
g = g[q]
ra = np.array(g['ra_cent'], dtype=np.float32)
dec = np.array(g['dec_cent'], dtype=np.float32)
gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs').galactic

cat = fits.open('starcat_allscans_06_12_19.fits')[1].data
w = np.where((cat['expsum'] > 4) & (cat['ctsum'] > 1))
cat = cat[w]
catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')
cid, gid, angsep, ang3d = gal.search_around_sky(catgal, 3*u.degree)
q = np.unique(cid)
carea = cat[q]


# Or load this instead
cat = fits.open('starcat_allscans_06_12_19.fits')[1].data
w = np.where((cat['expsum'] > 4) & (cat['ctsum'] > 1))
cat = cat[w]
carea = fits.open('plane_ingaisarea.fits')[1].data
r = np.where(carea['FLUXERR_AUTO']/carea['FLUX_AUTO'] < 0.2)
carea = carea[r]
galex = fits.open('GAISPlane.fits')[1].data
q = np.where((galex['nuv_magerr'] > 0) & (galex['nuv_magerr'] < 100) & (galex['nuv_magerr'] < 0.2) & (galex['primary_flag'] == 1) & (galex['nuv_weight'] > 0))
galex = galex[q]


plt.hist(galex['nuv_mag'], range=[12,25], bins=20, label='GAIS', log=True, alpha=0.7, color='#00E6FF')
plt.hist(cat['nuv'], bins=20, range=[12,25], histtype='step', linewidth=2, stacked=True, label='UVGAPS', log=True, color='black')
plt.hist(carea['nuv'], range=[12,25], bins=20, histtype='step', linewidth=2, stacked=True, label='UVGAPS in GAIS area', log=True, color='red')
plt.legend(scatterpoints=1, loc=1)
plt.xlabel('NUV')
plt.ylim(1000, 10**7)
plt.show()



############################################################
# How many objects in the plane survey are in none of the others = 483,956
############################################################
cat = fits.open('starcat_allscans_06_12_19.fits')[1].data # catalog
w = np.where((cat['expsum'] > 4) & (cat['ctsum'] > 1))
cat = Table(cat[w])
ps = fits.open('plane_ps2_09_01_19.fits')[1].data # Pan starrs
cg = fits.open('plane_gaiadr2_dust_06_12_19.fits')[1].data  # gaia
gg = fits.open('plane_gais_06_12_19.fits')[1].data # GAIS 


catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
gggal = SkyCoord(gg['ALPHA_J2000']*u.deg, gg['DELTA_J2000']*u.deg, frame='icrs')
cind, ggind, a, b = search_around_sky(catgal, gggal, 0.1*u.arcsec)
cat.remove_rows(cind)

catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
psgal = SkyCoord(ps['ALPHA_J2000']*u.deg, ps['DELTA_J2000']*u.deg, frame='icrs')
cind, pind, a, b = search_around_sky(catgal, psgal, 0.1*u.arcsec)
cat.remove_rows(cind)


catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
cggal = SkyCoord(cg['ra_plane']*u.deg, cg['dec_plane']*u.deg, frame='icrs')
cind, cgind, a, b = search_around_sky(catgal, cggal, 0.1*u.arcsec)
cat.remove_rows(cind)

len(cat)

############################################################
# Exposure time vs gb
############################################################
exp = fits.open('expdata_11_28_18.fits')[1].data

scans = ['scan0.5', 'scan1.4', 'scan2.3', 'scan3.2', 'scan4.1', 'scan5.0', 'scan5.9', 'scan6.8', 'scan8.6', 'scan9.5', 'scan10.4', 'scan11.3', 'scan12.2', 'scan14.0', 'scan14.9', 'scan15.8', 'scan16.7', 'scan17.6', 'scan18.5', 'scan19.4', 'scan20.3', 'scan21.2', 'scan22.1', 'scan23.0', 'scan23.9', 'scan24.8', 'scan25.7', 'scan28.4', 'scan29.3', 'scan30.2', 'scan31.1', 'scan32.0', 'scan32.9', 'scan33.8', 'scan34.7', 'scan35.6', 'scan39.2', 'scan42.8', 'scan43.7', 'scan44.6', 'scan45.5', 'scan46.4', 'scan47.3', 'scan48.2', 'scan49.1', 'scan50.0', 'scan67.1', 'scan68.9', 'scan71.6', 'scan74.3', 'scan75.2', 'scan76.1', 'scan77.0', 'scan77.9', 'scan78.8', 'scan79.7', 'scan80.6', 'scan81.5', 'scan82.4', 'scan83.3', 'scan87.8', 'scan88.7', 'scan89.6', 'scan90.5', 'scan91.4', 'scan92.3', 'scan93.2', 'scan94.1', 'scan95.0', 'scan95.9', 'scan96.8', 'scan97.7', 'scan98.6', 'scan99.5', 'scan100.4', 'scan101.3', 'scan102.2', 'scan103.1', 'scan104.0', 'scan104.9', 'scan105.8', 'scan106.7', 'scan107.6', 'scan110.3', 'scan111.2', 'scan112.1', 'scan113.0', 'scan113.9', 'scan114.8', 'scan119.3', 'scan121.1', 'scan122.9', 'scan124.7', 'scan125.6', 'scan126.5', 'scan127.4', 'scan128.3', 'scan129.2', 'scan130.1', 'scan131.0', 'scan131.9', 'scan132.8', 'scan133.7', 'scan134.6', 'scan135.5', 'scan136.4', 'scan137.3', 'scan138.2', 'scan139.1', 'scan140.0', 'scan140.9', 'scan141.8', 'scan143.6', 'scan144.5', 'scan145.4', 'scan148.1', 'scan149.0', 'scan149.9', 'scan150.8', 'scan151.7', 'scan152.6', 'scan153.5', 'scan155.3', 'scan156.2', 'scan157.1', 'scan158.0', 'scan160.7', 'scan161.6', 'scan163.4', 'scan167.0', 'scan167.9', 'scan172.4', 'scan173.3', 'scan174.2', 'scan175.1', 'scan176.0', 'scan176.9', 'scan177.8', 'scan178.7', 'scan179.6', 'scan180.5', 'scan183.2', 'scan185.0', 'scan190.4', 'scan191.3', 'scan197.6', 'scan198.5', 'scan200.3', 'scan201.2', 'scan203.0', 'scan203.9', 'scan205.7', 'scan206.6', 'scan207.5', 'scan208.4', 'scan209.3', 'scan210.2', 'scan211.1', 'scan212.0', 'scan212.9', 'scan213.8', 'scan214.7', 'scan215.6', 'scan216.5', 'scan217.4', 'scan218.3', 'scan219.2', 'scan220.1', 'scan221.0', 'scan221.9', 'scan222.8', 'scan223.7', 'scan224.6', 'scan225.5', 'scan226.4', 'scan228.2', 'scan229.1', 'scan230.0', 'scan230.9', 'scan231.8', 'scan234.5', 'scan235.4', 'scan236.3', 'scan237.2', 'scan238.1', 'scan239.0', 'scan239.9', 'scan240.8', 'scan241.7', 'scan242.6', 'scan243.5', 'scan244.4', 'scan245.3', 'scan246.2', 'scan247.1', 'scan248.0', 'scan248.9', 'scan249.8', 'scan250.7', 'scan251.6', 'scan252.5', 'scan253.4', 'scan254.3', 'scan255.2', 'scan256.1', 'scan257.0', 'scan258.8', 'scan259.7', 'scan260.6', 'scan261.5', 'scan263.3', 'scan264.2', 'scan265.1', 'scan266.0', 'scan266.9', 'scan268.7', 'scan269.6', 'scan270.5', 'scan271.4', 'scan272.3', 'scan273.2', 'scan274.1', 'scan275.0', 'scan275.9', 'scan276.8', 'scan278.6', 'scan279.5', 'scan281.3', 'scan283.1', 'scan284.0', 'scan285.8', 'scan286.7', 'scan288.5', 'scan289.4', 'scan290.3', 'scan291.2', 'scan292.1', 'scan293.0', 'scan293.9', 'scan295.7', 'scan297.5', 'scan298.4', 'scan301.1', 'scan302.0', 'scan302.9', 'scan303.8', 'scan304.7', 'scan305.6', 'scan306.5', 'scan308.3', 'scan309.2', 'scan310.1', 'scan315.5', 'scan316.4', 'scan317.3', 'scan318.2', 'scan319.1', 'scan320.0', 'scan320.9', 'scan321.8', 'scan322.7', 'scan323.6', 'scan324.5', 'scan325.4', 'scan326.3', 'scan327.2', 'scan328.1', 'scan329.0', 'scan329.9', 'scan331.7', 'scan332.6', 'scan333.5', 'scan334.4', 'scan335.3', 'scan338.0', 'scan338.9', 'scan339.8', 'scan341.6', 'scan342.5', 'scan343.4', 'scan345.2', 'scan348.8', 'scan349.7', 'scan350.6', 'scan351.5', 'scan352.4', 'scan353.3', 'scan354.2', 'scan355.1', 'scan356.0', 'scan357.8', 'scan358.7', 'scan359.0']

expvals = 0
for i in range(len(scans)):
    expvals += exp[scans[i]]
    #if (i % 20 == 0) & (i != 0):
    expvals = expvals#/20.
    plt.plot(np.linspace(-10, 10, len(expvals)), expvals, label=str(i), linewidth=0.1, alpha=0.3)
    expvals = 0
plt.xlabel('Galactic Latitude')
plt.ylabel('Exposure time [s]')
plt.ylim(0, 200)
plt.show()






a5 = a[:len(a)/4]
a6 = a[len(a)/4:len(a)/2]
a7 = a[len(a)/2:len(a)*3/4]
a8 = a[len(a)*3/4:]
plt.scatter(a['ra'], a['dec'], s=1)
plt.scatter(a5['ra'], a5['dec']-50, s=1, c='red')
plt.scatter(a6['ra'], a6['dec']-50, s=1, c='green')
plt.scatter(a7['ra'], a7['dec']-50, s=1, c='purple')
plt.scatter(a8['ra'], a8['dec']-50, s=1, c='orange')
plt.show()
ascii.write(a5, 'starcat_coords_pt5_07_03_19.txt', format='basic')
ascii.write(a6, 'starcat_coords_pt6_07_03_19.txt', format='basic')
ascii.write(a7, 'starcat_coords_pt7_07_03_19.txt', format='basic')
ascii.write(a8, 'starcat_coords_pt8_07_03_19.txt', format='basic')


addhash('starcat_coords_pt1_07_03_19.txt')
addhash('starcat_coords_pt2_07_03_19.txt')
addhash('starcat_coords_pt3_07_03_19.txt')
addhash('starcat_coords_pt4_07_03_19.txt')
addhash('starcat_coords_pt5_07_03_19.txt')
addhash('starcat_coords_pt6_07_03_19.txt')
addhash('starcat_coords_pt7_07_03_19.txt')
addhash('starcat_coords_pt8_07_03_19.txt')




ps = fits.open('plane_ps2_09_01_19.fits')[1].data
dra = (ps['ALPHA_J2000']-ps['raMean'])*3600
ddec = (ps['DELTA_J2000']-ps['decMean'])*3600
w = np.where((np.sqrt(dra**2+ddec**2) < 3) & (ps['ng'] > 5) & (ps['nr'] > 5) & (ps['ni'] > 5) & (ps['nz'] > 5) & (ps['ny'] > 5) & (ps['expsum'] > 4) & (ps['ctsum'] > 1))


