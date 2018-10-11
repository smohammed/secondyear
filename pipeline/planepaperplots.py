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
cat = fits.open('../starcat_allscans_09-19-18')[1].data # catalog
ps = fits.open('../plane_ps1_g10-20_09_26_18.fits')[1].data # Pan starrs
cg = fits.open('../plane_gaiadr2_dust_09_27_18.fits')[1].data  # gaia
gg = fits.open('../plane_gais_09_26_18.fits')[1].data # GAIS 

############################################################
# Now make plots
############################################################
# 1. Gaia CMD
gag = fits.open('gais_gaiadr2_negparcuts_dust_09_27_18.fits')[1].data 

negpar = np.where((cg['dist'] > 0) & (cg['visibility_periods_used'] > 8) & (cg['phot_bp_mean_mag'] > 0) & (cg['phot_rp_mean_mag'] > 0) & (cg['expsum'] > 10) & (cg['ebv'] > 0) & (cg['parallax_error']/cg['parallax'] < 0.1))
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

############################################################
# Angular separations for Gaia and PS1
############################################################
cg = fits.open('plane_gaiadr2_dust_09_27_18.fits')[1].data
plt.hist(cg['angsep']*3600, bins=20),plt.show()
plt.xlabel('Angular Separation [Plane to Gaia]')
plt.show()

ps = fits.open('plane_ps1_g10-20_09_26_18.fits')[1].data
plt.hist(ps['angsep']*3600, bins=20)
plt.xlabel('Angular Separation [Plane to PS1]')
plt.show()


############################################################
# NUV comparison
############################################################
gg = fits.open('plane_gais_09-26-18.fits')[1].data

averages = []
deln = gg['nuv_mag']-gg['nuv']
delnavg = gg['nuv_mag']-gg['nuv']
for mag in np.arange(12.5, 20.5, 0.5):
	cut = np.where((gg['nuv_mag'] > mag) & (gg['nuv_mag'] < mag+0.5))
	avg = np.average(deln[cut]) - 0.03
	delnavg[cut] = delnavg[cut] - avg
	averages.append(avg)

scatter_contour(gg['nuv_mag'], gg['nuv_mag']-gg['nuv'], threshold=10000, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))

plt.axhline(y=0, color='red')
plt.xlim(12, 22.5)
plt.ylim(-1, 1.5)
plt.xlabel('NUV (GAIS)')
plt.ylabel('$\Delta$NUV (GAIS - Plane)')
plt.show()




# For all data
averages = [-0.5088040321564006, -0.30568558997508744, -0.08878561654404307, 0.06906811684820328, 0.15215738430922907, 0.18585467134802813, 0.2044312305142438, 0.21806953210465735, 0.22379168603395064, 0.22521566946951102, 0.2222612442213893, 0.21768934311656635, 0.23189637728274717, 0.3016845341048976, 0.3764177765040423, 0.6407062596745744]
a = dict(zip(np.arange(12.5, 20.5, 0.5), averages))


# Modify NUV values
cat = fits.open('../starcat_allscans_09-19-18')[1].data # catalog

	for mag in np.arange(12.5, 20.5, 0.5):
		cut = np.where((matched['nuv'] > mag) & (matched['nuv'] < mag+0.5))
		matched['nuv'][cut] = matched['nuv'][cut] - a[mag]
	lower = np.where(matched['nuv'] < 12.5)
	higher = np.where(matched['nuv'] > 20.5)
	matched['nuv'][lower] = matched['nuv'][lower] - a[12.5]
	matched['nuv'][higher] = matched['nuv'][higher] - a[20]
#ascii.write(matched, '...', format='basic')

############################################################
# NUV histogram for all surveys
############################################################
cat = fits.open('starcat_allscans_09-28-18.fits')[1].data # catalog
ps = fits.open('plane_ps1_g10-20_09_26_18.fits')[1].data # Pan starrs
cg = fits.open('plane_gaiadr2_dust_09_27_18.fits')[1].data  # gaia
gg = fits.open('plane_gais_09-26-18.fits')[1].data # GAIS 
allcat = fits.open('plane_gaia_gais_ps1_all_09_27_18.fits')[1].data
none = Table.read('plane_innosurveys.txt', format='ascii')

'''
cat = cat[np.where(cat['nuv'] < 20.)]
ps = ps[np.where(ps['nuv'] < 20.)]
cg = cg[np.where(cg['nuv'] < 20.)]
gg = gg[np.where(gg['nuv'] < 20.)]
allcat = allcat[np.where(allcat['nuv'] < 20.)]
none = none[np.where(none['nuv'] < 20.)]
'''

ps = ps[np.where(ps['nuv'] < 20.5)]

plt.hist(cat['nuv'], bins=20, range=[12,21], label='Plane')
plt.hist(cg['nuv'], bins=20, range=[12, 21], histtype='step',linewidth=3, stacked=True, label='Plane+Gaia')
plt.hist(ps['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=3, stacked=True, label='Plane+PS1')
plt.hist(gg['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='Plane+GAIS')
plt.hist(none['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='No surveys', color='yellow')
plt.hist(allcat['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='Plane+All surveys')

plt.legend(loc=2)
plt.xlabel('NUV')
plt.show()

############################################################
# NUV offset histogram with GAIS
############################################################
for mag in range(12, 21,2):
    cut = np.where((gg['nuv'] > mag) & (gg['nuv'] < mag+2))
    dnuv = (gg['nuv_mag']-gg['nuv'])[cut] 
    plt.hist(dnuv, histtype='step', fill=False, stacked=True, label=str(mag)+'-'+str(mag+2), range=[-1, 1], bins=40)
plt.legend(scatterpoints=1, loc=2)
plt.xlabel('dNUV (GAIS - plane)')
plt.show()


############################################################
# g vs g-r
############################################################
threshold = 10000
bins = 4

scatter_contour(ps['g_ps']-ps['r_ps'], ps['g_ps'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=(bins)), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.gca().invert_yaxis()
plt.xlim(-7, 7.5)
plt.ylim(20, 11)
plt.xlabel('g - r')
plt.ylabel('g')
plt.show()


############################################################
# g-r vs NUV-g
############################################################
pickles = Table.read('picklemags_laphare_final.txt', format='ascii')
scatter_contour(ps['nuv']-ps['g_ps'], ps['g_ps']-ps['r_ps'], threshold=10000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.scatter(pickles['nuv']-pickles['g'], pickles['g']-pickles['r'], color='red', label='SED model',  s=30)
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5, 0.6), 2.5, 0.5, facecolor='yellow', alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5, 1.1), 6.5, 0.5, facecolor='yellow', alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((2, -0.1), 6, 0.8, facecolor='lightblue', alpha=0.5, angle=5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-1.9, -0.75), 3, 0.65, facecolor='gray', alpha=0.5))
plt.arrow(3,  -0.75,  2.972-1.1838,  1.1838-0.8664,  head_length=0.05,  head_width=0.02,  color='red')
plt.xlim((-2.1, 9.1))
plt.ylim((-1, 2))

plt.annotate('O', xy=(0, -0.9), size=20)
plt.annotate('B', xy=(1, -0.81), size=20)
plt.annotate('A', xy=(2.5, -0.5), size=20)
plt.annotate('F', xy=(3.7, -0.4), size=20)
plt.annotate('G', xy=(5, -0.2), size=20)
plt.annotate('K', xy=(6.2, -0.1), size=20)
plt.xlabel('NUV - g')
plt.ylabel('g - r')
plt.legend(scatterpoints=1, loc=4)
plt.show()

####################################################################
# NUV - g vs g - i
####################################################################
scatter_contour(ps['g_ps']-ps['i_ps'],ps['nuv']-ps['g_ps'],threshold=10000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.scatter(pickles['g']-pickles['i'],pickles['nuv']-pickles['g'],color='red', label='SED model', s=30, zorder=10)
plt.arrow(-1.9, 2, 1.1936-0.6533, 2.9720-1.1936, head_length=0.05, head_width=0.02, color='red', zorder=10)
plt.xlim((-2, 3))
plt.ylim((8, -2))

plt.annotate('O', xy=(-1.2, 0), size=20)
plt.annotate('B', xy=(-1.1, 1), size=20)
plt.annotate('A', xy=(-0.9, 2.3), size=20)
plt.annotate('F', xy=(-0.5, 3.5), size=20)
plt.annotate('G', xy=(1.1, 5), size=20)
plt.annotate('K', xy=(1.5, 6), size=20)

plt.xlabel('g - i')
plt.ylabel('NUV - g')
plt.legend(scatterpoints=1, loc=3)
plt.show()


############################################################
# Area
############################################################
g = fits.open('gr67tile_table.fits')[1].data
ra = np.array(g['ra_cent'], dtype=np.float32)
dec = np.array(g['dec_cent'], dtype=np.float32)
gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
gl = gal.galactic.l.degree
gb = gal.galactic.b.degree
q = np.where((gb > -10) & (gb < 10))
g = g[q]
ra = np.array(g['ra_cent'], dtype=np.float32)
dec = np.array(g['dec_cent'], dtype=np.float32)
gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs').galactic

cat = fits.open('starcat_allscans_09-28-18.fits')[1].data
catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')
cid, gid, angsep, ang3d = gal.search_around_sky(catgal, 1*u.degree)
q = np.unique(cid)
c2 = cat[q]

gais = fits.open('GALEXPlanenew.fits')[1].data

plt.hist(gais['nuv_mag'], range=[12,22], bins=20, label='GAIS', log=True)
plt.hist(c2['nuv'], range=[12,22], bins=20, histtype='step', linewidth=2, stacked=True, label='Plane in GAIS area', log=True)
plt.legend(scatterpoints=1, loc=2)
plt.xlabel('NUV')
plt.show()



############################################################
# How many objects in the plane survey are in none of the others = 483,956
############################################################
cat = fits.open('starcat_allscans_09-19-18.fits')[1].d28a # catalog
ps = fits.open('plane_ps1_g10-20_09_26_18.fits')[1].data # Pan starrs
cg = fits.open('plane_gaiadr2_dust_09_27_18.fits')[1].data  # gaia
gg = fits.open('plane_gais_09-26-18.fits')[1].data # GAIS 


catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
psgal = SkyCoord(ps['ALPHA_J2000']*u.deg, ps['DELTA_J2000']*u.deg, frame='icrs')
cggal = SkyCoord(cg['ALPHA_J2000']*u.deg, cg['DELTA_J2000']*u.deg, frame='icrs')
gggal = SkyCoord(gg['ALPHA_J2000']*u.deg, gg['DELTA_J2000']*u.deg, frame='icrs')


cat = Table(cat)

cind, ggind, a, b = search_around_sky(catgal, gggal, 0.1*u.arcsec)
cat.remove_rows(cind)

catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
psgal = SkyCoord(ps['ALPHA_J2000']*u.deg, ps['DELTA_J2000']*u.deg, frame='icrs')
cind, pind, a, b = search_around_sky(catgal, psgal, 0.1*u.arcsec)
cat.remove_rows(cind)


catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
cind, cgind, a, b = search_around_sky(catgal, cggal, 0.1*u.arcsec)
cat.remove_rows(cind)

len(cat)
