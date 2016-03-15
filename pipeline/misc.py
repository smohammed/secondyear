####################################################################
# Rename 2MASS catalog columns
####################################################################
star.rename_column('cntr_01','cntr')
star.rename_column('number_01','number')
star.rename_column('x_image_01','x_image')
star.rename_column('y_image_01','y_image')
star.rename_column('flux_auto_01','flux_auto')
star.rename_column('fluxerr_auto_01','fluxerr_auto')
star.rename_column('x_new_01','x_new')
star.rename_column('y_new_01','y_new')
star.rename_column('nuv_01','nuv')
star.rename_column('gl_01','gl_sex')
star.rename_column('gb_01','gb_sex')
star.rename_column('ra_01','ra_sex')
star.rename_column('dec_01','dec_sex')
star.rename_column('ra','ra_2mass')
star.rename_column('dec','dec_2mass')
star.rename_column('j_m','j')
star.rename_column('h_m','h')
star.rename_column('k_m','k')

####################################################################
# Import different Jlim and arcsec search radii from 2MASS
####################################################################
a2125 = Table.read('newfield_gal_ipac_2mass_matches_2arcsec_j_lt12.5_tests.txt',format='ascii')
a213 = Table.read('newfield_gal_ipac_2mass_matches_2arcsec_j_lt13_tests.txt',format='ascii')
a2135 = Table.read('newfield_gal_ipac_2mass_matches_2arcsec_j_lt13.5_tests.txt',format='ascii')
a214 = Table.read('newfield_gal_ipac_2mass_matches_2arcsec_j_lt14_tests.txt',format='ascii')
a2145 = Table.read('newfield_gal_ipac_2mass_matches_2arcsec_j_lt14.5_tests.txt',format='ascii')
b3125 = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt12.5_tests.txt',format='ascii')
b313 = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt13_tests.txt',format='ascii')
b3135 = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt13.5_tests.txt',format='ascii')
b314 = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt14_tests.txt',format='ascii')
b3145 = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt14.5_tests.txt',format='ascii')
c4125 = Table.read('newfield_gal_ipac_2mass_matches_4arcsec_j_lt12.5_tests.txt',format='ascii')
c413 = Table.read('newfield_gal_ipac_2mass_matches_4arcsec_j_lt13_tests.txt',format='ascii')
c4135 = Table.read('newfield_gal_ipac_2mass_matches_4arcsec_j_lt13.5_tests.txt',format='ascii')
c414 = Table.read('newfield_gal_ipac_2mass_matches_4arcsec_j_lt14_tests.txt',format='ascii')
c4145 = Table.read('newfield_gal_ipac_2mass_matches_4arcsec_j_lt14.5_tests.txt',format='ascii')


#Plot radius search vs J mag search limits 
files = [a2125,a213,a2135,a214,a2145,b3125,b313,b3135,b314,b3145,c4125,c413,c4135,c414,c4145]
jlim = [12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5]
arcsec = [2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]

# Field density plot vs J mag search limits
field = (20*7.4)*3600.
den, denlo,denhi = [],[],[]
for i in files:
    den.append(len(i)/field)    
    denlo.append(len(i[np.where(np.abs(i['gb_sex']) < 5.)])/field)
    denhi.append(len(i[np.where(np.abs(i['gb_sex']) > 5.)])/field)

# For whole field
plt.scatter(jlim,den,c=arcsec,edgecolor='none',s=40)
plt.xlabel('Jlim [mag]')
plt.ylabel('Counts/area [#/arcsec$^2$]')
cm = plt.colorbar()
cm.set_label('Limiting radius [arcsec]')
for i in range(len(files)):
    plt.annotate(str(len(files[i])),xy=(jlim[i]+0.05,den[i]))
plt.show()

# Now broken up by lower/upper galactic plane
plt.scatter(jlim,denlo,c=arcsec,edgecolor='none',s=80,marker='s',label='abs(b) < 5')
plt.scatter(jlim,denhi,c=arcsec,edgecolor='none',s=80,marker='o',label='abs(b) > 5')
plt.xlabel('Jlim [mag]')
plt.ylabel('Counts/area [#/arcsec$^2$]')
cm = plt.colorbar()
cm.set_label('Limiting radius [arcsec]')7
for i in range(len(files)):
    plt.annotate(str(len(files[i][np.where(np.abs(files[i]['gb_sex']) < 5.)])),xy=(jlim[i]+0.05,denlo[i]))
    plt.annotate(str(len(files[i][np.where(np.abs(files[i]['gb_sex']) > 5.)])),xy=(jlim[i]+0.05,denhi[i]))

plt.title('2MASS match counts by J_lim and radius search')
plt.legend(scatterpoints=1,loc=2)
plt.show()

####################################################################
# Add WCS coords to plots
####################################################################
img = fits.open('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
galex = Table.read('galex0data.txt',format='ascii')
gallim = np.where((galex['gl_galex'] > 0) & (galex['gl_galex'] < 7.5) & (galex['gb_galex']> -10) & (galex['gb_galex'] < 10))
galex = galex[gallim]

from astropy.wcs import WCS
header = fits.getheader('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')
w = WCS(header)
x0,y0 = w.wcs_pix2world(0,0,0)
xn,yn = w.wcs_pix2world(4860,12840,0)

xn = xn - 360

matplotlib.rcParams['figure.figsize'] = 7.5, 20
plt.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray,extent=[x0,xn,y0,yn])

#plt.scatter(star['gl_sex'],star['gb_sex'],edgecolor='red',facecolor='none',s=20)
plt.scatter(galex['gl_galex'],galex['gb_galex'],edgecolor='red',facecolor='red',s=5)
plt.scatter(vphas['l'],vphas['b'],edgecolor='blue',facecolor='blue',s=5,alpha=0.1)
plt.gca().invert_xaxis()
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.xlim((0,7.4))
plt.ylim((-10,10))
plt.show()

plt.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray)

plt.scatter(star['x_new'],star['y_new'],edgecolor='red',facecolor='none',s=20)

plt.show()

####################################################################
# Sky plot
####################################################################
from astropy.wcs import WCS
img = fits.open('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
galex = Table.read('galex0data.txt',format='ascii')
gallim = np.where((galex['gl_galex'] > 0) & (galex['gl_galex'] < 7.5) & (galex['gb_galex']> -10) & (galex['gb_galex'] < 10))
galex = galex[gallim]
header = fits.getheader('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')
w = WCS(header)
x0,y0 = w.wcs_pix2world(0,0,0)
xn,yn = w.wcs_pix2world(4860,12840,0)
xn = xn - 360

vphas = Table.read('vphas_gl_0_to_7.fits',format='fits')
vpgal = SkyCoord(vphas['RAJ2000']*u.deg,vphas['DEJ2000']*u.deg,frame='icrs')
galexgal = SkyCoord(galex['ra']*u.deg,galex['dec']*u.deg,frame='icrs')
galind, vpasind, angsep, dist3d = search_around_sky(galexgal,vpgal,2*u.arcsec)
galex2 = galex[galind]
vgalex = vphas[vpasind]
vgalex = hstack([galex2,vgalex])

plt.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray,extent=[x0,xn,y0,yn])
plt.scatter(galex['gl_galex'],galex['gb_galex'],edgecolor='red',facecolor='red',s=5,label='GALEX')
plt.scatter(vphas['l'],vphas['b'],edgecolor='blue',facecolor='blue',s=5,label='VPHAS')
plt.scatter(vgalex['l'],vgalex['b'],edgecolor='yellow',facecolor='yellow',s=5,label='VPHAS+GALEX')
plt.legend(loc=3,scatterpoints=1,prop={'size':12})
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.xlim((0,7.45))
plt.ylim((-10,10))
plt.gca().invert_xaxis()
plt.show()

from astropy.wcs import WCS
img = fits.open('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')[0].data
galex = Table.read('galex0data.txt',format='ascii')
gallim = np.where((galex['gl_galex'] > 0) & (galex['gl_galex'] < 7.5) & (galex['gb_galex']> -10) & (galex['gb_galex'] < 10))
galex = galex[gallim]
header = fits.getheader('newfield/count_map_05-68_gPr_cata_10_corr_gal.fits')
w = WCS(header)
x0,y0 = w.wcs_pix2world(0,0,0)
xn,yn = w.wcs_pix2world(4860,12840,0)
xn = xn - 360

vgalex = Table.read('galex0vphas.txt',format='ascii')

f, ((ax1, ax2)) = plt.subplots(1, 2, sharex='col')
ax1.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray,extent=[x0,xn,y0,yn])
ax1.set_xlabel('Galactic Longitude')
ax1.set_ylabel('Galactic Latitude')
ax1.invert_xaxis()
ax1.set_xlim((0,7.45))
ax1.set_ylim((-10,10))
ax2.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray,extent=[x0,xn,y0,yn])
ax2.scatter(galex['gl_galex'],galex['gb_galex'],edgecolor='red',facecolor='red',s=5,label='GALEX')
ax2.scatter(vphas['l'],vphas['b'],edgecolor='blue',facecolor='blue',s=5,label='VPHAS')
ax2.scatter(vgalex['l'],vgalex['b'],edgecolor='yellow',facecolor='yellow',s=5,label='VPHAS+GALEX')
ax2.set_xlabel('Galactic Longitude')
ax2.set_xlim((0,7.45))
ax2.set_ylim((-10,10))
#ax2.invert_xaxis()
ax2.get_yaxis().set_visible(False)
f.subplots_adjust(wspace=0)
plt.show()

######################################################################
# Reddening values 
######################################################################
NUV = 2.9720
B = 1.3429
V = 1.0
J = 0.2876
H = 0.1783
K = 0.1170
u = 1.5916
g = 1.1838
r = 0.8664
i = 0.6418

####################################################################
# Match GAIS + VPHAS 
####################################################################
galex = fits.open('GALEXAIS.fits')[1].data
vphas = fits.open('vphas_allg_gl0-40.fits')[1].data
vpcut = np.where(vphas['u_AB']-vphas['g_AB'] < 1.25)
vphas = vphas[vpcut]
galexgal = SkyCoord(galex['ra']*u.deg,galex['dec']*u.deg,frame='icrs')
vphasgal = SkyCoord(vphas['RAJ2000']*u.deg,vphas['DEJ2000']*u.deg,frame='icrs')
galexind,vphasind, angsep,sep3d = search_around_sky(galexgal,vphasgal,3.5*u.arcsec)

g2 = Table(galex[galexind])
v2 = Table(vphas[vphasind])
tot = hstack([g2,v2])
tot['angsep'] = angsep

print 'part 1 done'
vp1 = fits.open('vphas_allg_gl200-250.fits')[1].data
vp1cut = np.where(vp1['u_AB']-vp1['g_AB'] < 1.25)
vp1 = vp1[vp1cut]
vpgal = SkyCoord(vp1['RAJ2000']*u.deg,vp1['DEJ2000']*u.deg,frame='icrs')
galexind,vpind, angsep1,sep3d = search_around_sky(galexgal,vpgal,3.5*u.arcsec)
g2 = Table(galex[galexind])
v12 = Table(vp1[vpind])
tot1 = hstack([g2,v12])
tot1['angsep'] = angsep1

print 'part 2 done'
vp1 = fits.open('vphas_allg_gl250-300.fits')[1].data
vp1cut = np.where(vp1['u_AB']-vp1['g_AB'] < 1.25)
vp1 = vp1[vp1cut]
vpgal = SkyCoord(vp1['RAJ2000']*u.deg,vp1['DEJ2000']*u.deg,frame='icrs')
galexind,vpind, angsep2,sep3d = search_around_sky(galexgal,vpgal,3.5*u.arcsec)
g2 = Table(galex[galexind])
v12 = Table(vp1[vpind])
tot2 = hstack([g2,v12])
tot2['angsep'] = angsep2

print 'part 3 done'
vp1 = fits.open('vphas_allg_gl300-360.fits')[1].data
vp1cut = np.where(vp1['u_AB']-vp1['g_AB'] < 1.25)
vp1 = vp1[vp1cut]
vpgal = SkyCoord(vp1['RAJ2000']*u.deg,vp1['DEJ2000']*u.deg,frame='icrs')
galexind,vpind, angsep3,sep3d = search_around_sky(galexgal,vpgal,3.5*u.arcsec)
g2 = Table(galex[galexind])
v12 = Table(vp1[vpind])
tot3 = hstack([g2,v12])
tot3['angsep'] = angsep3
totall = vstack([tot,tot1,tot2,tot3])

####################################################################
# SED plot
####################################################################
gv = Table.read('GAIS_VPHAS.txt',format='ascii')

nuvi = wd2m['nuv_mag'] - wd2m['i_AB']
ui = wd2m['u_AB'] - wd2m['i_AB']
gi = wd2m['g_AB'] - wd2m['i_AB']
ri = wd2m['r_AB'] - wd2m['i_AB']
ji = wd2m['j'] - wd2m['i_AB']
hi = wd2m['h'] - wd2m['i_AB']
ki = wd2m['k'] - wd2m['i_AB']

cut = np.where(wd2m['j'] > 0)
nuvi2 = nuvi[cut]
ui2 = ui[cut]
gi2 = gi[cut]
ri2 = ri[cut]
ji2 = ji[cut]
hi2 = hi[cut]
ki2 = ki[cut]

cols = np.ones(len(wd2m)*7)
cols[0+len(wd2m)*0:len(wd2m)*1] = nuvi
cols[0+len(wd2m)*1:len(wd2m)*2] = ui
cols[0+len(wd2m)*2:len(wd2m)*3] = gi
cols[0+len(wd2m)*3:len(wd2m)*4] = ri
cols[0+len(wd2m)*4:len(wd2m)*5] = ji
cols[0+len(wd2m)*5:len(wd2m)*6] = hi
cols[0+len(wd2m)*6:len(wd2m)*7] = ki
cols2 = np.ones(len(wd2m[cut])*7)
cols2[0+len(wd2m[cut])*0:len(wd2m[cut])*1] = nuvi2
cols2[0+len(wd2m[cut])*1:len(wd2m[cut])*2] = ui2
cols2[0+len(wd2m[cut])*2:len(wd2m[cut])*3] = gi2
cols2[0+len(wd2m[cut])*3:len(wd2m[cut])*4] = ri2
cols2[0+len(wd2m[cut])*4:len(wd2m[cut])*5] = ji2
cols2[0+len(wd2m[cut])*5:len(wd2m[cut])*6] = hi2
cols2[0+len(wd2m[cut])*6:len(wd2m[cut])*7] = ki2


labels = ['NUV-i', 'u-i', 'g-i', 'r-i', 'J-i', 'H-i', 'K-i']
plt.xticks([1, 2, 3, 4, 5, 6, 7], labels)
for i in range(len(wd2m)):
    plt.plot([1, 2, 3, 4, 5, 6, 7], cols[[i, i+len(wd2m), i+len(wd2m)*2, i+len(wd2m)*3, i+len(wd2m)*4, i+len(wd2m)*5, i+len(wd2m)*6]], alpha=0.01, color='black')
for i in range(len(wd2m[cut])):
    plt.plot([1, 2, 3, 4, 5, 6, 7], cols2[[i, i+len(wd2m[cut]), i+len(wd2m[cut])*2, i+len(wd2m[cut])*3, i+len(wd2m[cut])*4, i+len(wd2m[cut])*5, i+len(wd2m[cut])*6]], alpha=0.05, color='red')

plt.plot([1, 2, 3, 4, 5, 6, 7], cols[[1, 1+len(wd2m), 1+len(wd2m)*2, 1+len(wd2m)*3, 1+len(wd2m)*4, 1+len(wd2m)*5, 1+len(wd2m)*6]], alpha=0.1, color='black', label='VPHAS')
plt.plot([1, 2, 3, 4, 5, 6, 7], cols2[[1, 1+len(wd2m[cut]), 1+len(wd2m[cut])*2, 1+len(wd2m[cut])*3, 1+len(wd2m[cut])*4, 1+len(wd2m[cut])*5, 1+len(wd2m[cut])*6]], alpha=0.1, color='red', label='VPHAS + 2MASS')
plt.legend(loc=2)
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda$ - i (ABmag)')
plt.ylim((5, -12))
plt.title('VPHAS + 2MASS,  Verbeek WDs')
plt.show()

####################################################################
# GALEX and VPHAS sky plot
####################################################################
galex = fits.open('GALEXAIS.fits')[1].data
galexcut = np.where((galex['glon'] > 0) & (galex['glon'] < 40))
g2 = galex[galexcut]
galexgal = SkyCoord(g2['glon']*u.deg, g2['glat']*u.deg, frame='galactic')
vp = fits.open('vphas_allg_gl0-40.fits')[1].data
vpgal = SkyCoord(vp['RAJ2000']*u.deg,vp['DEJ2000']*u.deg,frame='icrs').galactic
galexind,vpind,angsep,ang3d = search_around_sky(galexgal,vpgal,3*u.arcsec)

g3 = g2[galexind]

plt.scatter(g2['glon'],g2['glat'],edgecolor='none',label='galex')
plt.scatter(vpgal.l.degree,vpgal.b.degree,edgecolor='none',facecolor='red',label='vphas')
plt.scatter(g3['glon'],g3['glat'],edgecolor='none',facecolor='yellow',label='both')
plt.legend(scatterpoints=1)
plt.show()

galexcut = np.where((galex['glon'] > 200) & (galex['glon'] < 250))
g2 = galex[galexcut]
galexgal = SkyCoord(g2['glon']*u.deg,g2['glat']*u.deg,frame='galactic')
vp = fits.open('vphas_allg_gl200-250.fits')[1].data
vpgal = SkyCoord(vp['RAJ2000']*u.deg,vp['DEJ2000']*u.deg,frame='icrs').galactic
galexind,vpind,angsep,ang3d = search_around_sky(galexgal,vpgal,3*u.arcsec)
g3 = g2[galexind]
plt.scatter(g2['glon'],g2['glat'],edgecolor='none',label='galex')
plt.scatter(vpgal.l.degree,vpgal.b.degree,edgecolor='none',facecolor='red',label='vphas')
plt.scatter(g3['glon'],g3['glat'],edgecolor='none',facecolor='yellow',label='both')
plt.legend(scatterpoints=1)
plt.show()

galexcut = np.where((galex['glon'] > 250) & (galex['glon'] < 300))
g2 = galex[galexcut]
galexgal = SkyCoord(g2['glon']*u.deg,g2['glat']*u.deg,frame='galactic')
vp = fits.open('vphas_allg_gl250-300.fits')[1].data
vpgal = SkyCoord(vp['RAJ2000']*u.deg,vp['DEJ2000']*u.deg,frame='icrs').galactic
galexind,vpind,angsep,ang3d = search_around_sky(galexgal,vpgal,3*u.arcsec)
g3 = g2[galexind]
plt.scatter(g2['glon'],g2['glat'],edgecolor='none',label='galex')
plt.scatter(vpgal.l.degree,vpgal.b.degree,edgecolor='none',facecolor='red',label='vphas')
plt.scatter(g3['glon'],g3['glat'],edgecolor='none',facecolor='yellow',label='both')
plt.legend(scatterpoints=1)
plt.show()

galexcut = np.where((galex['glon'] > 300) & (galex['glon'] < 360))
g2 = galex[galexcut]
galexgal = SkyCoord(g2['glon']*u.deg,g2['glat']*u.deg,frame='galactic')
vp = fits.open('vphas_allg_gl300-360.fits')[1].data
vpgal = SkyCoord(vp['RAJ2000']*u.deg,vp['DEJ2000']*u.deg,frame='icrs').galactic
galexind,vpind,angsep,ang3d = search_around_sky(galexgal,vpgal,3*u.arcsec)
g3 = g2[galexind]
plt.scatter(g2['glon'],g2['glat'],edgecolor='none',label='galex')
plt.scatter(vpgal.l.degree,vpgal.b.degree,edgecolor='none',facecolor='red',label='vphas')
plt.scatter(g3['glon'],g3['glat'],edgecolor='none',facecolor='yellow',label='both')
plt.legend(scatterpoints=1)
plt.show()

fig,axes = plt.subplots(2,2,sharey=True)
xco = gv['glon']
yco = gv['g_AB'] - gv['r_AB']
axes[0,0].scatter(xco,yco,edgecolor='none')
axes[0,1].scatter(xco,yco,edgecolor='none')
axes[1,1].scatter(xco,yco,edgecolor='none')
axes[1,0].scatter(xco,yco,edgecolor='none')

wdx = wds['gl_galex']
wdy = wds['g_AB'] - wds['r_AB']
axes[0,0].scatter(wdx,wdy,c='red',label='All WDs')
axes[0,1].scatter(wdx,wdy,c='red')
axes[1,1].scatter(wdx,wdy,c='red')
axes[1,0].scatter(wdx,wdy,c='red')
'''
axes[0,0].scatter(gl[cut],gb[cut],c='green',label='NUV-g > 3')
axes[0,1].scatter(gl[cut],gb[cut],c='green')
axes[1,1].scatter(gl[cut],gb[cut],c='green')
axes[1,0].scatter(gl[cut],gb[cut],c='green')
'''

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
axes[1,0].set_ylabel('g - r')
axes[0,0].set_ylabel('g - r')
axes[0,0].legend(scatterpoints=1)
fig.subplots_adjust(wspace=0)
plt.show()

####################################################################
# g - r vs NUV - g with boxes
####################################################################
scatter_contour(gv['nuv_mag']-gv['g_AB'],gv['g_AB']-gv['r_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,0.6),2.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,1.1),6.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((2,0.1),6,0.6,facecolor='red',alpha=0.5,angle=5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-2,-0.75),3,0.65,facecolor='gray',alpha=0.5))
plt.xlim((-3,9))
plt.ylim((-1,3))
plt.show()

####################################################################
# g vs g-r by gl
####################################################################
v1 = fits.open('vphas_allg_gl0-40.fits')[1].data
m = 8./1.5
b = 22 - m*0.5
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 0 < gl < 40')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_0-40.png')
plt.clf()

v1cut = np.where(v1['u_AB']-v1['g_AB'] < 1.)
v1 = v1[v1cut]
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 0 < gl < 40, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_0-40_ug1cut.png')
plt.clf()

v1 = fits.open('vphas_allg_gl200-250.fits')[1].data
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 200 < gl < 250')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_200-250.png')
plt.clf()

v1cut = np.where(v1['u_AB']-v1['g_AB'] < 1.)
v1 = v1[v1cut]
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 200 < gl < 250, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_200-250_ug1cut.png')
plt.clf()

v1 = fits.open('vphas_allg_gl250-300.fits')[1].data
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 250 < gl < 300')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_250-300.png')
plt.clf()
v1cut = np.where(v1['u_AB']-v1['g_AB'] < 1.)
v1 = v1[v1cut]
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 250 < gl < 300, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_250-300_ug1cut.png')
plt.clf()

v1 = fits.open('vphas_allg_gl300-360.fits')[1].data
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 300 < gl < 360')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_300-360.png')
plt.clf()
v1cut = np.where(v1['u_AB']-v1['g_AB'] < 1.)
v1 = v1[v1cut]
scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 300 < gl < 360, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_300-360_ug1cut.png')
plt.clf()

######################################################################
# lambda - i SEDs
######################################################################
v1 = Table.read('wds_vphasonly.txt',format='ascii')
m = 8/1.5
b = 22-8/1.5*0.5
#cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))

v1cut = np.where((v1['gl'] > 300) & (v1['gl'] < 360))
v1 = v1[v1cut]

wd = Table.read('wds_gais_vphas_newgrcut.txt',format='ascii')
wdcut = np.where((wd['gl_galex'] > 300) & (wd['gl_galex'] < 360))
wd2 = wd[wdcut]

pwds = Table.read('picklemags_wds.txt',format='ascii')

# SED plot
v1nuvi = np.zeros(len(v1))
v1ui = v1['u_AB'] - v1['r_AB']
v1gi = v1['g_AB'] - v1['r_AB']
#v1ri = v1['r_AB'] - v1['r_AB']

nuvi = wd2['nuv_mag'] - wd2['r_AB']
ui = wd2['u_AB'] - wd2['r_AB']
gi = wd2['g_AB'] - wd2['r_AB']
#ri = wd2['r_AB'] - wd2['r_AB']

colsv1 = np.ones(len(v1)*3)
colsv1[0+len(v1)*0:len(v1)*1] = v1nuvi
colsv1[0+len(v1)*1:len(v1)*2] = v1ui
colsv1[0+len(v1)*2:len(v1)*3] = v1gi
#colsv1[0+len(v1)*3:len(v1)*4] = v1ri

cols = np.ones(len(wd2)*3)
cols[0+len(wd2)*0:len(wd2)*1] = nuvi
cols[0+len(wd2)*1:len(wd2)*2] = ui
cols[0+len(wd2)*2:len(wd2)*3] = gi
#cols[0+len(wd2)*3:len(wd2)*4] = ri

# All VPHAS

for i in range(len(v1)):
    plt.plot([1,2,3],colsv1[[i,i+len(v1),i+len(v1)*2]],alpha=0.05,color='blue')

# VPHAS + GAIS
for i in range(len(wd2)):
    plt.plot([1,2,3],cols[[i,i+len(wd2),i+len(wd2)*2]],alpha=0.1,color='red')

# Plot this to make legend
plt.plot([1,2,3],colsv1[[1,1+len(v1),1+len(v1)*2]],alpha=0.4,color='black',label='VPHAS only')
plt.plot([1,2,3],cols[[1,1+len(wd2),1+len(wd2)*2]],alpha=0.4,color='red',label='VPHAS + GAIS')

plt.plot([1,2,3],[pwds['nuv']-pwds['r'],pwds['u']-pwds['r'],pwds['g']-pwds['r']],color='black')
plt.plot([1],[pwds['nuv'][0]-pwds['r'][0]],color='black',label='SED')
plt.legend(loc=4)
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda$ - i (ABmag)')
plt.ylim((4,-2))
plt.title('VPHAS + GAIS WDs, new gr cut, 300 < gl < 360')
labels = ['NUV-r','u-r','g-r']
plt.xticks([1,2,3],labels)

plt.show()

# Or, Pickles SED points with only GAIS+VPHAS points
for i in range(len(pickles)):
    plt.plot([1,2,3,4],[pickles[i]['nuv']-pickles[i]['i'],pickles[i]['u']-pickles[i]['i'],pickles[i]['g']-pickles[i]['i'],pickles[i]['r']-pickles[i]['i']],color='blue',alpha=0.3)
# VPHAS + GAIS
for i in range(len(wd2)):
    plt.plot([1,2,3,4],cols[[i,i+len(wd2),i+len(wd2)*2,i+len(wd2)*3]],alpha=0.1,color='red')

# Plot this to make legend
plt.plot([1,2,3,4],cols[[1,1+len(wd2),1+len(wd2)*2,1+len(wd2)*3]],alpha=0.4,color='red',label='VPHAS + GAIS')
plt.plot([1,2,3,4],[pickles[1]['nuv']-pickles[1]['i'],pickles[1]['u']-pickles[1]['i'],pickles[1]['g']-pickles[1]['i'],pickles[1]['r']-pickles[1]['i']],color='blue',alpha=0.4,label='Pickles')
plt.legend(loc=3)
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda$ - i (ABmag)')
plt.ylim((3,-1.5))
plt.title('VPHAS + GAIS WDs, new gr cut, 300 < gl < 360')
labels = ['NUV-i','u-i','g-i','r-i']
plt.xticks([1,2,3,4],labels)
plt.show()


######################################################################
# lambda - i histograms
######################################################################
v = Table.read('wds_vphasonly.txt',format='ascii')
wd = Table.read('wds_gais_vphas_newgrcut.txt',format='ascii')
pwds = Table.read('picklemags_wds.txt',format='ascii')

v1cut = np.where((v['gl'] > 0) & (v['gl'] < 40))
v2cut = np.where((v['gl'] > 200) & (v['gl'] < 250))
v3cut = np.where((v['gl'] > 250) & (v['gl'] < 300))
v4cut = np.where((v['gl'] > 300) & (v['gl'] < 360))
wd1cut = np.where((wd['gl_galex'] > 0) & (wd['gl_galex'] < 40))
wd2cut = np.where((wd['gl_galex'] > 200) & (wd['gl_galex'] < 250))
wd3cut = np.where((wd['gl_galex'] > 250) & (wd['gl_galex'] < 300))
wd4cut = np.where((wd['gl_galex'] > 300) & (wd['gl_galex'] < 360))

bins= np.linspace(-2,3,20)
fig,axes = plt.subplots(2,2,sharey=True,sharex=True)
axes[0,0].hist((v['r_AB']-v['i_AB'])[v1cut],bins=bins,label='VPHAS')
axes[0,1].hist((v['r_AB']-v['i_AB'])[v2cut],bins=bins)
axes[1,0].hist((v['r_AB']-v['i_AB'])[v3cut],bins=bins)
axes[1,1].hist((v['r_AB']-v['i_AB'])[v4cut],bins=bins)
axes[0,0].hist((wd['r_AB']-wd['i_AB'])[wd1cut],bins=bins, label='VPHAS+GAIS')
axes[0,1].hist((wd['r_AB']-wd['i_AB'])[wd2cut],bins=bins)
axes[1,0].hist((wd['r_AB']-wd['i_AB'])[wd3cut],bins=bins)
axes[1,1].hist((wd['r_AB']-wd['i_AB'])[wd4cut],bins=bins)
axes[0,0].axvline(x=pwds['r'][0]-pwds['i'][0],c='black',linewidth=2,label='SED')
axes[0,1].axvline(x=pwds['r'][0]-pwds['i'][0],c='black',linewidth=2)
axes[1,0].axvline(x=pwds['r'][0]-pwds['i'][0],c='black',linewidth=2)
axes[1,1].axvline(x=pwds['r'][0]-pwds['i'][0],c='black',linewidth=2)
axes[0,0].annotate('0 < gl < 40',xy=(1,400),size=15)
axes[0,1].annotate('200 < gl < 250',xy=(1,400),size=15)
axes[1,0].annotate('250 < gl < 300',xy=(1,400),size=15)
axes[1,1].annotate('300 < gl < 360',xy=(1,400),size=15)
fig.subplots_adjust(wspace=0,hspace=0)
fig.suptitle(r'WDCs')
axes[1,0].set_xlabel('r - i (ABmag)')
axes[1,1].set_xlabel('r - i (ABmag)')
axes[0,0].legend(scatterpoints=1)
plt.show()

# NUV only
bins= np.linspace(-3,5,20)
fig,axes = plt.subplots(2,2,sharey=True,sharex=True)
axes[0,0].hist((wd['nuv_mag']-wd['i_AB'])[wd1cut],bins=bins, label='VPHAS+GAIS',color='green')
axes[0,1].hist((wd['nuv_mag']-wd['i_AB'])[wd2cut],bins=bins,color='green')
axes[1,0].hist((wd['nuv_mag']-wd['i_AB'])[wd3cut],bins=bins,color='green')
axes[1,1].hist((wd['nuv_mag']-wd['i_AB'])[wd4cut],bins=bins,color='green')
axes[0,0].annotate('0 < gl < 40',xy=(2,80),size=15)
axes[0,1].annotate('200 < gl < 250',xy=(2,80),size=15)
axes[1,0].annotate('250 < gl < 300',xy=(2,80),size=15)
axes[1,1].annotate('300 < gl < 360',xy=(2,80),size=15)
axes[0,0].axvline(x=pwds['nuv'][0]-pwds['i'][0],c='black',linewidth=2,label='avg SED')
axes[0,1].axvline(x=pwds['nuv'][0]-pwds['i'][0],c='black',linewidth=2)
axes[1,0].axvline(x=pwds['nuv'][0]-pwds['i'][0],c='black',linewidth=2)
axes[1,1].axvline(x=pwds['nuv'][0]-pwds['i'][0],c='black',linewidth=2)

fig.subplots_adjust(wspace=0,hspace=0)

#fig.suptitle(r'WDs, $\Delta \theta_{max}$ = 1.5", gr cut')
fig.suptitle(r'WDCs')
axes[1,0].set_xlabel('NUV - i (ABmag)')
axes[1,1].set_xlabel('NUV - i (ABmag)')
axes[0,0].legend(scatterpoints=1)
plt.show()


######################################################################
# lambda - i histograms WITH DUST
######################################################################
vphas = Table.read('wds_vphasonly.txt', format='ascii')
wd = Table.read('wds_gais_vphas_newgrcut.txt', format='ascii')
vpcut = np.where((vphas['gl'] > 200) & (vphas['gl'] < 250))
wdcut = np.where((wd['gl_galex'] > 200) & (wd['gl_galex'] < 250))
vphas = vphas[vpcut][:5000]
wd = wd[wdcut]


pwds = Table.read('picklemags_wds.txt', format='ascii')
# VPHAS + WD
nuv_wdav = Table.read('wd_GV_av1_dm_1ER_nuv.txt', format='ascii')
u_wdav = Table.read('wd_GV_av1_dm_1ER_u.txt', format='ascii')
g_wdav = Table.read('wd_GV_av1_dm_1ER_g.txt', format='ascii')
r_wdav = Table.read('wd_GV_av1_dm_1ER_r.txt', format='ascii')
i_wdav = Table.read('wd_GV_av1_dm_1ER_i.txt', format='ascii')
# VPHAS only
u_vpav = Table.read('wd_V_av1_dm_1ER_u.txt', format='ascii')
g_vpav = Table.read('wd_V_av1_dm_1ER_g.txt', format='ascii')
r_vpav = Table.read('wd_V_av1_dm_1ER_r.txt', format='ascii')
i_vpav = Table.read('wd_V_av1_dm_1ER_i.txt', format='ascii')

temp = '10k'

fig,axes = plt.subplots(2,2,sharex=True)
bins= np.linspace(-2,3,20)
#axes[0,0].hist((vphas['_AB']-vphas['i_AB'])[v1cut],bins=bins)
axes[0,1].hist((vphas['u_AB']-u_vpav[temp])-(vphas['i_AB']-i_vpav[temp]),bins=bins,label='VP')
axes[1,0].hist((vphas['g_AB']-g_vpav[temp])-(vphas['i_AB']-i_vpav[temp]),bins=bins)
axes[1,1].hist((vphas['r_AB']-r_vpav[temp])-(vphas['i_AB']-i_vpav[temp]),bins=bins)

axes[0,0].hist((wd['nuv_mag']-nuv_wdav[temp])-(wd['i_AB']-i_wdav[temp]),bins=bins)
axes[0,1].hist((wd['u_AB']-u_wdav[temp])-(wd['i_AB']-i_wdav[temp]),bins=bins, label='VP+G')
axes[1,0].hist((wd['g_AB']-g_wdav[temp])-(wd['i_AB']-i_wdav[temp]),bins=bins)
axes[1,1].hist((wd['r_AB']-r_wdav[temp])-(wd['i_AB']-i_wdav[temp]),bins=bins)


axes[0,0].axvline(x=pwds['nuv'][0]-pwds['i'][0],c='black',linewidth=2)
axes[0,1].axvline(x=pwds['u'][0]-pwds['i'][0],c='black',linewidth=2,label='SED')
axes[1,0].axvline(x=pwds['g'][0]-pwds['i'][0],c='black',linewidth=2)
axes[1,1].axvline(x=pwds['r'][0]-pwds['i'][0],c='black',linewidth=2)
axes[0,0].annotate('NUV - i',xy=(1,80),size=15)
axes[0,1].annotate('u - i',xy=(1,300),size=15)
axes[1,0].annotate('g - i',xy=(1,400),size=15)
axes[1,1].annotate('r - i',xy=(1,400),size=15)
fig.subplots_adjust(hspace=0)
fig.suptitle(r'WDCs, 200 < gl < 250, with dust corr, 10k K')

axes[0,0].set_ylim((0,180))
axes[0,1].set_ylim((0,1400))
axes[1,0].set_ylim((0,1600))
axes[1,1].set_ylim((0,2500))

axes[1,0].set_xlabel('$\lambda$ - i (ABmag)')
axes[1,1].set_xlabel('$\lambda$ - i (ABmag)')
axes[0,1].legend(scatterpoints=1)
plt.show()



plt.scatter(vphas['i_AB'],i_vpav['10k'],label='10k')
plt.scatter(vphas['i_AB'],i_vpav['15k'],label='15k',c='red')
plt.scatter(vphas['i_AB'],i_vpav['20k'],label='20k',c='green')
plt.scatter(vphas['i_AB'],i_vpav['25k'],label='25k',c='orange',alpha=0.5)
plt.scatter(vphas['i_AB'],i_vpav['30k'],label='30k',c='purple',alpha=0.5)
plt.legend(scatterpoints=1,loc=2)
plt.xlabel('i (AB)')
plt.ylabel('A$_{V}$')
plt.title('WDC, VP+G')
plt.ylim((0,4))
plt.show()


temp = '20k'
scatter_contour(vp['g_AB']-vp['r_AB'],vp['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
plt.scatter((vphas['g_AB']-g_vpav[temp])-(vphas['r_AB']-r_vpav[temp]),vphas['g_AB']-g_vpav[temp],alpha=0.3)
plt.xlim((-1,3))
plt.ylim((23,12))
plt.xlabel('g - r (AB)')
plt.ylabel('g (AB)')
plt.title('VPHAS + WDs, 200 < gl < 250, dust, 20k K')
plt.show()


v1cut = np.where(vphas['u_AB']-vphas['g_AB'] < 1.)
v1 = v1[v1cut]

scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
cut = np.where(((v1['g_AB']-v1['r_AB'])*m+b < v1['g_AB']) & (v1['g_AB'] > 19.))
plt.scatter((v1['g_AB']-v1['r_AB'])[cut],v1['g_AB'][cut],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((22,12))
plt.title('VPHAS only, 0 < gl < 40, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
plt.savefig('11-25-gvsgr_vphasonly_0-40_ug1cut.png')
plt.clf()

v1cut = np.where(v1['u_AB']-v1['g_AB'] < 1.)

cut = np.where(( v1['g_AB'] > (v1['g_AB']-v1['r_AB'])*m+b) & (v1['g_AB'] > 19.))


cut = np.where((vpgal.galactic.l.degree > 0.) & (vpgal.galactic.l.degree < 40.))
v1 = v1[cut]

wdcut = np.where((vpwdgal.l.degree > 0.) & (vpwdgal.l.degree < 40.))
wdv1 = vpwd[wdcut]

scatter_contour((v1['g_AB']-v1['r_AB']),v1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))

plt.scatter((wdv1['g_AB']-wdv1['r_AB']),wdv1['g_AB'],edgecolor='none',facecolor='red')
plt.xlim((-1,3))
plt.ylim((23,12))
plt.title('VPHAS only, 0 < gl < 40')#, u - g < 1')
plt.xlabel('g - r')
plt.ylabel('g')
#plt.savefig('11-25-gvsgr_vphasonly_0-40_ug1cut.png')
#plt.clf()


galex = fits.open('GALEXAIS.fits')[1].data
vpgal = SkyCoord(vp['RAJ2000']*u.deg,vp['DEJ2000']*u.deg,frame='icrs')
galexgal = SkyCoord(galex['ra']*u.deg,galex['dec']*u.deg,frame='icrs')
galexind,vpind,angsep,sep3d = search_around_sky(galexgal,vpgal,3.5*u.arcsec)

gv1 = np.where((gv['gl_galex'] > 0.) & (gv['gl_galex'] < 40.))
gv2 = np.where((gv['gl_galex'] > 200.) & (gv['gl_galex'] < 250.))
gv3 = np.where((gv['gl_galex'] > 250.) & (gv['gl_galex'] < 300.))
gv4 = np.where((gv['gl_galex'] > 300.) & (gv['gl_galex'] < 360.))
gv1 = gv[gv1]
gv2 = gv[gv2]
gv3 = gv[gv3]
gv4 = gv[gv4]


gv1ug = np.where(gv1['u_AB']-gv1['g_AB'] < 1.)
gv2ug = np.where(gv2['u_AB']-gv2['g_AB'] < 1.)
gv3ug = np.where(gv3['u_AB']-gv3['g_AB'] < 1.)
gv4ug = np.where(gv4['u_AB']-gv4['g_AB'] < 1.)

gv1ug = gv1[gv1ug]
gv2ug = gv2[gv2ug]
gv3ug = gv3[gv3ug]
gv4ug = gv4[gv4ug]

scatter_contour(gv3ug['g_AB']-gv3ug['r_AB'],gv3ug['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),contour_args=dict(),plot_args=dict(color='k',markersize=1))
plt.xlabel('g - r (ABmag)')
plt.ylabel('g (ABmag)')
plt.title('GAIS+VPHAS, 300 < gl < 360, u - g < 1')
plt.xlim((-1,3))
plt.ylim((23,12))
plt.show()


skyrange = ['1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '260.6', '261.5', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '289.4', '290.3', '291.2', '292.1', '293.0', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

skyrange = ['102.2', '103.1', '104.0', '104.9', '10.4', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '11.3', '114.8', '119.3', '12.2', '14.0', '14.9', '1.4', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '2.3', '24.8', '25.7', '271.4', '272.3', '273.2', '274.1', '278.6', '279.5', '281.3', '283.1', '284.0', '28.4', '289.4', '290.3', '291.2', '292.1', '293.0', '29.3', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '30.2', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '31.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '32.0', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '32.9', '3.2', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '33.8', '339.8', '341.6', '342.5', '343.4', '345.2', '34.7', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '35.6', '357.8', '358.7', '359.6', '39.2', '4.1', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '5.0', '5.9', '5', '67.1', '68.9', '6.8', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '8.6', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '9.5', '96.8', '97.7', '98.6', '99.5']
data = Table.read('starcat_5.txt',format='ascii')

for region in skyrange:
    t1 = Table.read('starcat_'+region+'.txt',format='ascii')
    data = vstack([data, t1])
    print region

ascii.write(data,'starcat_1-4.txt',format='basic')

skyrange = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '260.6', '261.5', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '289.4', '290.3', '291.2', '292.1', '293.0', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

count = []
for region in skyrange:
    t1 = Table.read('starcat_'+region+'.txt', format='ascii')
    count.append(len(t1))

comb1.rename_column('nuv','nuv_sex')
comb1.rename_column('gl','gl_sex')
comb1.rename_column('gb','gb_sex')
comb1.rename_column('ra','ra_sex')
comb1.rename_column('dec','dec_sex')
comb1.rename_column('Glon','gl_t2')
comb1.rename_column('Glat','gb_t2')
comb1.rename_column('RAJ2000','ra_t2')
comb1.rename_column('DEJ2000','dec_t2')

cat1 = fits.open('Dunmaps/starcat_total_1-5.fits')[1].data
cat1gal = SkyCoord(cat1['gl']*u.deg,cat1['gb']*u.deg,frame='galactic')

cat2 = fits.open('Dunmaps/starcat_total_1-4.fits')[1].data
cat2gal = SkyCoord(cat2['gl']*u.deg,cat2['gb']*u.deg,frame='galactic')

tycho = fits.open('tycho2.fits')[1].data
tychogal = SkyCoord(tycho['Glon']*u.deg,tycho['Glat']*u.deg,frame='galactic')

for region in range(0,360,5):
    cat1 = comb1[np.where((comb1['gl_t2'] > region) & (comb1['gl_t2'] < region+5))]
    dgl1 = (cat1['gl_sex']-cat1['gl_t2'])
    dgb1 = (cat1['gb_sex']-cat1['gb_t2'])

    cat2 = comb2[np.where((comb2['gl_t2'] > region) & (comb2['gl_t2'] < region+5))]
    dgl2 = (cat2['gl_sex']-cat2['gl_t2'])
    dgb2 = (cat2['gb_sex']-cat2['gb_t2'])

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.scatter(dgl1*3600, dgb1*3600, edgecolor='none', alpha=0.1)
    ax1.axhline(y=0,c='black')
    ax1.axvline(x=0,c='black')
    ax2.scatter(dgl2*3600, dgb2*3600, edgecolor='none', alpha=0.1)
    ax2.axhline(y=0,c='black')
    ax2.axvline(x=0,c='black')

    ax1.set_xlabel('$\Delta$ gl')
    ax2.set_xlabel('$\Delta$ gl')
    ax1.set_ylabel('$\Delta$ gb')
    ax1.set_title('old, gl = '+str(region)+'-'+str(region+5)+', avg=('+str(np.mean(dgl1)*3600)[:4]+', '+str(np.mean(dgb1)*3600)[:4]+')')
    ax2.set_title('new, gl = '+str(region)+'-'+str(region+5)+', avg=('+str(np.mean(dgl2)*3600)[:4]+', '+str(np.mean(dgb2)*3600)[:4]+')')
    fig.subplots_adjust(wspace=0)

    plt.savefig('03-11-sex-t2_offset_'+str(region)+'-'+str(region+5)+'.png')
    plt.clf()
    print region

####################################################################
#check sextractor nuv filter comparisons
####################################################################

sky = ['4.1', '5.0', '19.4', '18.5', '328.1', '329.0']
galex = fits.open('../../GALEXAIS.fits')[1].data
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')

convlist = ['gauss2.0_5x5', 'gauss2.5_5x5', 'gauss3.0_5x5', 'gauss1.5_3x3', 'gauss2.0_3x3', 'gauss3.0_7x7', 'gauss4.0_7x7', 'gauss5.0_9x9', 'mexhat1.5_5x5', 'mexhat2.0_7x7', 'mexhat2.5_7x7', 'mexhat3.0_9x9', 'mexhat4.0_9x9', 'mexhat5.0_11x11', 'tophat1.5_3x3', 'tophat2.0_3x3', 'tophat2.5_3x3', 'tophat3.0_3x3', 'tophat4.0_5x5', 'tophat5.0_5x5']

for conv in convlist:
    for region in sky:
        a20 = Table.read('starcat_'+region+'_'+conv+'.txt', format='ascii')
        norm = Table.read('../starcat_'+region+'.txt', format='ascii')

        a20gal = SkyCoord(a20['gl']*u.deg, a20['gb']*u.deg, frame='galactic')
        normgal = SkyCoord(norm['gl']*u.deg, norm['gb']*u.deg, frame='galactic')

        a20ind,  galex20ind,  angsep20,  ang3d = search_around_sky(a20gal, galexgal, 3.5*u.arcsec)
        normind,  galexnormind,  angsepnorm,  ang3d = search_around_sky(normgal, galexgal, 3.5*u.arcsec)

        a20 = a20[a20ind]
        norm = norm[normind]
        g20 = Table(galex[galex20ind])
        gnorm = Table(galex[galexnormind])

        comb20 = hstack([a20, g20])
        combnorm = hstack([norm, gnorm])

        comb20['angsep'] = angsep20
        combnorm['angsep'] = angsepnorm
        ascii.write(comb20, 'sextractor_galex_matches'+region+'_'+conv+'.txt', format='basic')

        fig, (ax0, ax1) = plt.subplots(2, sharex=True)
        ax0.scatter(combnorm['nuv_mag'], combnorm['nuv']-combnorm['nuv_mag'], edgecolor='none', alpha=0.1)
        ax0.axhline(y=0, c='black')
        ax1.scatter(comb20['nuv_mag'], comb20['nuv']-comb20['nuv_mag'], edgecolor='none', alpha=0.1)
        ax1.axhline(y=0, c='black')
        ax1.set_xlabel('NUV$_{GAIS}$')
        ax0.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
        ax0.set_title('gl ='+str(region)+', default')
        ax0.set_title('gl ='+str(region)+', '+conv)
        ax0.set_xlim((12, 23))
        ax0.set_ylim((-4, 2))
        ax1.set_xlim((12, 23))
        ax1.set_ylim((-4, 2))
        ax0.annotate('Default, N = '+str(len(combnorm)), xy=(13.5, -3))
        ax1.annotate(conv+', N = '+str(len(comb20)), xy=(13.5, -3))
        fig.subplots_adjust(hspace=0)
        plt.savefig('03-13-nuvcomp_gl'+region+conv+'.png')
        plt.clf()

####################################################################
# Find duplicates in scans
####################################################################

a = Table.read('starcat_1.4.txt', format='ascii')
b = Table.read('starcat_2.3.txt', format='ascii')
agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
bgal = SkyCoord(b['gl']*u.deg, b['gb']*u.deg, frame='galactic')
aind, bind, angsep, ang3d = search_around_sky(agal, bgal, 1*u.arcsec)
a2 = a[aind]
b2 = b[bind]

plt.scatter(a['gl'], a['gb'], edgecolor='none')
plt.scatter(b['gl'], b['gb'], c='red', edgecolor='none')
plt.scatter(a2['gl'], a2['gb'], color='green', alpha=0.5)
plt.scatter(b2['gl'], b2['gb'], color='orange', alpha=0.5)

