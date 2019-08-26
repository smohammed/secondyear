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
plt.xlim(12, 21)
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
cg = fits.open('plane_gaiadr2_dust_06_12_19.fits')[1].data
ps = fits.open('plane_ps2_09_01_19.fits')[1].data
ps = Table(ps)
dra = (ps['ALPHA_J2000']-ps['raMean'])*3600
ddec = (ps['DELTA_J2000']-ps['decMean'])*3600
angsep = np.sqrt(dra**2+ddec**2)
ps['angsep'] = angsep

cg = cg[np.where(cg['angsep']*3600. < 2.0)]
ps = ps[np.where(ps['angsep'] < 2.0)]

plt.hist(cg['angsep']*3600, bins=20, alpha=0.7, label='Gaia match', color='#00E6FF', edgecolor='black')
plt.hist(ps['angsep'], bins=20, alpha=0.7, label='PS2 match', color='#FF5C5C', edgecolor='black')
plt.xlabel('Angular Separation')
plt.legend(scatterpoints=1)
plt.show()


############################################################
# NUV comparison
############################################################
gg = fits.open('plane_gais_11_27_18.fits')[1].data

'''
averages = []
deln = gg['nuv_mag']-gg['nuv_plane']
delnavg = gg['nuv_mag']-gg['nuv_plane']
for mag in np.arange(12.5, 20.5, 0.5):
	cut = np.where((gg['nuv_mag'] > mag) & (gg['nuv_mag'] < mag+0.5))
	avg = np.average(deln[cut]) - 0.03
	delnavg[cut] = delnavg[cut] - avg
	averages.append(avg)
'''

# Subtract 0.25 mag from plane NUV because the offset was added to the _cut.fits file
scatter_contour(gg['nuv_mag'], gg['nuv_mag']-(gg['nuv']-0.25), threshold=10000, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))

plt.axhline(y=0, color='red')
plt.xlim(12, 22.5)
plt.ylim(-1, 1.5)
plt.xlabel('NUV (GAIS)')
plt.ylabel('$\Delta$NUV (GAIS - Plane)')
plt.show()


############################################################
# NUV histogram for all surveys
############################################################
cat = fits.open('starcat_allscans_10-12-18_cuts.fits')[1].data # catalog
ps = fits.open('plane_ps2_09_01_19.fits')[1].data # Pan starrs
cg = fits.open('plane_gaiadr2_dust_11_14_18.fits')[1].data  # gaia
gg = fits.open('plane_gais_11_27_18.fits')[1].data # GAIS 
allcat = fits.open('plane_gaiadr2_gais_ps1_11_27_18.fits')[1].data
none = fits.open('plane_innosurveys_11_27_18.fits')[1].data
#cg = cg[np.where(cg['angsep']*3600. < 2.0)]
#ps = ps[np.where((ps['angsep']*3600. < 2.0) & (ps['expsum'] > 3.) & (ps['ctsum'] > 1))]
#ps = ps[~np.isnan(cat['bkgdsum'])]


plt.hist(cat['nuv'], bins=20, range=[12,21], label='UVGAPS', alpha=0.7, color='#00E6FF')
plt.hist(none['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='No surveys', color='yellow')
plt.hist(cg['nuv'], bins=20, range=[12, 21], histtype='step',linewidth=3, stacked=True, label='UVGAPS+Gaia')
plt.hist(gg['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='UVGAPS+GAIS')
plt.hist(ps['nuv'], bins=20, range=[12, 21], histtype='step', linewidth=3, stacked=True, label='UVGAPS+PS1')
plt.hist(allcat['nuv_plane'], bins=20, range=[12, 21], histtype='step', linewidth=2, stacked=True, label='UVGAPS+All surveys')

plt.legend(loc=2)
plt.xlabel('NUV')
plt.show()

############################################################
# NUV offset histogram with GAIS
############################################################
gg = fits.open('plane_gais_11_27_18.fits')[1].data # GAIS 

nuv = gg['nuv'] - 0.25

for mag in range(12, 21,2):
    cut = np.where((nuv > mag) & (nuv < mag+2))
    dnuv = (gg['nuv_mag']-nuv)[cut] 
    plt.hist(dnuv, histtype='step', fill=False, stacked=True, label=str(mag)+'-'+str(mag+2), range=[-1, 1], bins=40)
plt.legend(scatterpoints=1, loc=2)
plt.xlabel('$\Delta$NUV (GAIS - plane)')
plt.xlim(-0.4, 1)
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

#cat = fits.open('starcat_allscans_10-12-18_cuts.fits')[1].data
catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')
cid, gid, angsep, ang3d = gal.search_around_sky(catgal, 3*u.degree)
q = np.unique(cid)
carea = cat[q]

#galex = fits.open('GAISPlane.fits')[1].data

plt.hist(galex['nuv_mag'], range=[12,25], bins=20, label='GAIS', log=True, alpha=0.7, color='#00E6FF')
plt.hist(carea['nuv'], range=[12,25], bins=20, histtype='step', linewidth=2, stacked=True, label='UVGAPS in GAIS area', log=True, color='red')
plt.legend(scatterpoints=1, loc=2)
plt.xlabel('NUV')
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


