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
scatter_contour(vcat['nuv_sex']-vcat['g_AB'],vcat['g_AB']-vcat['r_AB'],threshold=1100,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
plt.scatter(vwd['nuv_sex']-vwd['g_AB'],vwd['g_AB']-vwd['r_AB'], edgecolor='none', facecolor='blue', label='WDCs')
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen', label='pickles', s=30)
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,0.6),2.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,1.1),6.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((2,0.1),6,0.6,facecolor='lightblue',alpha=0.5,angle=5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-2,-0.75),3,0.65,facecolor='gray',alpha=0.5))
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.xlim((-4,9))
plt.ylim((-1,3))

plt.annotate('O', xy=(0, -0.9), size=20)
plt.annotate('B', xy=(1, -0.81), size=20)
plt.annotate('A', xy=(2.5, -0.5), size=20)
plt.annotate('F', xy=(3.7, -0.4), size=20)
plt.annotate('G', xy=(5, -0.2), size=20)
plt.annotate('K', xy=(6.2, -0.1), size=20)

plt.xlabel('NUV - g (ABmag)')
plt.ylabel('g - r (ABmag)')
plt.legend(scatterpoints=1, loc=2)
plt.title('vPHAS + SExtractor')
plt.show()

####################################################################
# NUV - g vs g - i
####################################################################
scatter_contour(vcat['g_AB']-vcat['i_AB'],vcat['nuv_sex']-vcat['g_AB'],threshold=1100,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
plt.scatter(vwd['g_AB']-vwd['i_AB'],vwd['nuv_sex']-vwd['g_AB'], edgecolor='none', facecolor='blue', label='WDCs')
plt.scatter(pickles['g']-pickles['i'],pickles['nuv']-pickles['g'],color='darkgreen', label='pickles', s=30)
plt.arrow(-1.9, 2, 1.1936-0.6533, 2.9720-1.1936, head_length=0.05, head_width=0.02, color='red')
plt.xlim((-2, 3))
plt.ylim((8, -2))

plt.annotate('O', xy=(-1.2, 0), size=20)
plt.annotate('B', xy=(-1.1, 1), size=20)
plt.annotate('A', xy=(-0.9, 2.3), size=20)
plt.annotate('F', xy=(-0.4, 4.5), size=20)
plt.annotate('G', xy=(-0.25, 6), size=20)
plt.annotate('K', xy=(-0.1, 7), size=20)

plt.xlabel('g - i (ABmag)')
plt.ylabel('NUV - g (ABmag)')
plt.legend(scatterpoints=1, loc=3)
plt.title('vPHAS + SExtractor')
plt.show()


####################################################################
# g - r vs u - g
####################################################################
scatter_contour(vcat['u_AB']-vcat['g_AB'],vcat['g_AB']-vcat['r_AB'],threshold=1500,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
plt.scatter(vwd['u_AB']-vwd['g_AB'],vwd['g_AB']-vwd['r_AB'], edgecolor='none', facecolor='blue', label='WDCs')
plt.scatter(pickles['u']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen', label='pickles', s=30)
plt.arrow(3, -0.5, 1.5812 - 1.1936, 1.1936 - 0.8694, head_length=0.05, head_width=0.02, color='red')
plt.xlim((-1, 5))
plt.ylim((-1, 3))

plt.annotate('O', xy=(-0.3, -0.8), size=20)
plt.annotate('B', xy=(0.25, -0.7), size=20)
plt.annotate('A', xy=(1.05, -0.5), size=20)
plt.annotate('F', xy=(1.3, -0.55), size=20)
plt.annotate('G', xy=(1.8, -0.5), size=20)
plt.annotate('K', xy=(2.3, 0.0), size=20)
plt.annotate('M', xy=(2.8, 0.4), size=20)

plt.xlabel('u - g (ABmag)')
plt.ylabel('g - r (ABmag)')
plt.legend(scatterpoints=1, loc=2)
plt.title('vPHAS + SExtractor')
plt.show()


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
# Label all scans
####################################################################
s1 = ['1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4']

'101.3', 
s2 = ['110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5']

'200.3', 

s3 = ['201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']

table = Table.read('starcat_200.3mapweight.txt', format='ascii')
n = list(itertools.repeat('200.3', len(table)))
table['region'] = n

for region in s3:
    print region
    a = Table.read('starcat_'+region+'mapweight.txt',format='ascii')
    n = list(itertools.repeat(region, len(a)))
    a['region'] = n
    table = vstack((table, a))

ascii.write(table, 'starcat_set200-360.txt',format='basic')


galex = fits.open('../GALEXAIS.fits')[1].data
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')
filt = 'mexhat_4.0_9x9'
code = 'det_thresh2.5_phot_autopar2.5_3.5'

fits.HDUList([fits.PrimaryHDU(table)]).writeto('sex_galex_total.fits')


for region in t:
    cat = Table.read('Dunmaps/'+region,format='ascii')
    sexgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')
    sexind,  galexind,  angsep,  ang3d = search_around_sky(sexgal, galexgal, 3.5*u.arcsec)
    comb = ([cat[sexind],Table(galex[galexind])])
    comb['angsep'] = angsep
    comb.rename_column('nuv','nuv_sex')
    comb.rename_column('nuv_mag','nuv_galex')
    comb.rename_column('gl','gl_sex')
    comb.rename_column('gb','gb_sex')
    comb.rename_column('glon','gl_galex')
    comb.rename_column('glat','gb_galex')
    comb.rename_column('ra_1','ra_sex')
    comb.rename_column('dec_1','dec_sex')
    comb.rename_column('ra_2','ra_galex')
    comb.rename_column('dec_2','dec_galex')

    ascii.write(comb,'sex_galex_matches_'+region+'.txt',format='basic')



cat = Table.read('sex_tycho_matches_total_mapweight.txt', format='ascii')
c1 = cat[np.where(cat['gb_tycho'] < 0)]
c2 = cat[np.where(cat['gb_tycho'] > 0)]

for region in range(0,360,5):
    cat1 = c1[np.where((c1['gl_tycho'] > region) & (c1['gl_tycho'] < region+5))]
    dgl1 = (cat1['gl_sex']-cat1['gl_tycho'])
    dgb1 = (cat1['gb_sex']-cat1['gb_tycho'])

    cat2 = c2[np.where((c2['gl_tycho'] > region) & (c2['gl_tycho'] < region+5))]
    dgl2 = (cat2['gl_sex']-cat2['gl_tycho'])
    dgb2 = (cat2['gb_sex']-cat2['gb_tycho'])

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.scatter(dgl1*3600, dgb1*3600, edgecolor='none', alpha=0.1)
    ax1.axhline(y=0,c='black')
    ax1.axvline(x=0,c='black')
    ax2.scatter(dgl2*3600, dgb2*3600, edgecolor='none', alpha=0.1)
    ax2.axhline(y=0,c='black')
    ax2.axvline(x=0,c='black')

    ax1.set_xlabel('gl$_{SEx}$ - gl$_{Tycho2}$ [arcsec]')
    ax2.set_xlabel('gl$_{SEx}$ - gl$_{Tycho2}$ [arcsec]')
    ax1.set_ylabel('gb$_{SEx}$ - gb$_{Tycho2}$ [arcsec]')
    ax1.set_xlim((-4,4))
    ax1.set_ylim((-4,4))
    ax2.set_xlim((-4,4))
    ax2.set_ylim((-4,4))
    ax1.annotate('N$_{lower}$ = '+str(len(cat1)), xy=(-3.9,-3.9))
    ax2.annotate('N$_{upper}$ = '+str(len(cat2)), xy=(-3.9,-3.9))
    ax1.set_title('gb<0, gl='+str(region)+'-'+str(region+5)+', avg=('+str(np.mean(dgl1)*3600)[:4]+', '+str(np.mean(dgb1)*3600)[:4]+')')
    ax2.set_title('gb>0, avg=('+str(np.mean(dgl2)*3600)[:4]+', '+str(np.mean(dgb2)*3600)[:4]+')')
    fig.subplots_adjust(wspace=0)

    plt.savefig('04-13-sex-tycho_offset_'+str(region)+'-'+str(region+5)+'.png')
    plt.clf()
    print region


agal = SkyCoord(a['gl']*u.deg,a['gb']*u.deg, frame='galactic')
aind,  gaind,  angsepa,  ang3d = search_around_sky(agal, galexgal, 3.5*u.arcsec)
a2 = Table(a[aind])
g2 = Table(galex[gaind])
comb = hstack([a2, g2])
comb['angsep'] = angsepa
comb.rename_column('nuv', 'nuv_sex')
comb.rename_column('nuv_mag', 'nuv_galex')
comb.rename_column('gl', 'gl_sex')
comb.rename_column('gb', 'gb_sex')
comb.rename_column('ra_1', 'ra_sex')
comb.rename_column('dec_1', 'dec_sex')
comb.rename_column('ra_2', 'ra_galex')
comb.rename_column('dec_2', 'dec_galex')
comb.rename_column('glon', 'gl_galex')
comb.rename_column('glat', 'gb_galex')
combcut = np.where(comb['nuv_galex'] == -999.)
comb.remove_rows(combcut)

plt.scatter(comb['nuv_galex'], comb['nuv_sex']-comb['nuv_galex'], alpha=0.1, edgecolor='none')
plt.axhline(y=0, c='green')
plt.xlabel('NUV$_{GAIS}$')
plt.ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
plt.title('gl = 329.0')
plt.xlim((12, 23))
plt.ylim((-3, 2))
plt.xlim((12, 23))
plt.ylim((-3, 2))
plt.annotate('N = '+str(len(comb)), xy=(13, -2.5))
plt.show()

for x in range(0,100):
    for y in range(0,20):
        cut = np.where((cat['gl'] > x) & (cat['gl'] < x+1) & (cat['gb'] > y-10) & (cat['gb'] < y-10+1))
        cmap[x][y] = len(cat[cut])



cat = Table.read('sex_vphas_matches_total_mapweight.txt', format='ascii')



for region in range(0,360,10):
    cat1 = c1[np.where((c1['gl_tycho'] > region) & (c1['gl_tycho'] < region+5))]
    dgl1 = (cat1['gl_sex']-cat1['gl_tycho'])
    dgb1 = (cat1['gb_sex']-cat1['gb_tycho'])

    cat2 = c2[np.where((c2['gl_tycho'] > region) & (c2['gl_tycho'] < region+5))]
    dgl2 = (cat2['gl_sex']-cat2['gl_tycho'])
    dgb2 = (cat2['gb_sex']-cat2['gb_tycho'])

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.scatter(dgl1*3600, dgb1*3600, edgecolor='none', alpha=0.1)
    ax1.axhline(y=0,c='black')
    ax1.axvline(x=0,c='black')
    ax2.scatter(dgl2*3600, dgb2*3600, edgecolor='none', alpha=0.1)
    ax2.axhline(y=0,c='black')
    ax2.axvline(x=0,c='black')

    ax1.set_xlabel('gl$_{SEx}$ - gl$_{Tycho2}$ [arcsec]')
    ax2.set_xlabel('gl$_{SEx}$ - gl$_{Tycho2}$ [arcsec]')
    ax1.set_ylabel('gb$_{SEx}$ - gb$_{Tycho2}$ [arcsec]')
    ax1.set_xlim((-4,4))
    ax1.set_ylim((-4,4))
    ax2.set_xlim((-4,4))
    ax2.set_ylim((-4,4))
    ax1.annotate('N$_{lower}$ = '+str(len(cat1)), xy=(-3.9,-3.9))
    ax2.annotate('N$_{upper}$ = '+str(len(cat2)), xy=(-3.9,-3.9))
    ax1.set_title('gb<0, gl='+str(region)+'-'+str(region+5)+', avg=('+str(np.mean(dgl1)*3600)[:4]+', '+str(np.mean(dgb1)*3600)[:4]+')')
    ax2.set_title('gb>0, avg=('+str(np.mean(dgl2)*3600)[:4]+', '+str(np.mean(dgb2)*3600)[:4]+')')
    fig.subplots_adjust(wspace=0)

    plt.savefig('04-13-sex-tycho_offset_'+str(region)+'-'+str(region+5)+'.png')
    plt.clf()
    print region


for i in range(0, 360, 60):
    c2 = cat[np.where((cat['nuv'] < 19) & (cat['gl'] > i) & (cat['gl'] < i+60))]
    plt.hist2d(c2['gl'],c2['gb'], bins=[600,240], vmin=0,vmax=20, cmap=cm.gray)
    plt.title('SExtractor, NUV < 19, 600x240')
    plt.xlim((i, i+60))
    plt.ylim((-10,10))
    plt.savefig('04-18-sky_'+str(i)+'-'+str(i+60)+'sex.png')
    plt.clf()
    
    cg2 = gcat[np.where((gcat['nuv_sex'] < 19) & (gcat['nuv_galex'] < 19) & (gcat['gl_sex'] > i) & (gcat['gl_sex'] < i+60))]
    plt.hist2d(cg2['gl_sex'],cg2['gb_sex'], bins=[600,240], vmin=0,vmax=20, cmap=cm.gray)
    plt.title('SEx+GALEX, NUV < 19, 600x240')
    plt.xlim((i, i+60))
    plt.ylim((-10,10))
    plt.savefig('04-18-sky_'+str(i)+'-'+str(i+60)+'galex_sex.png')
    plt.clf()

    gal2 = galex[np.where((galex['nuv_mag'] < 19) & (galex['glon'] > i) & (galex['glon'] < i+60))]
    plt.hist2d(gal2['glon'],gal2['glat'], bins=[600,240], vmin=0,vmax=20, cmap=cm.gray)
    plt.title('GALEX, NUV < 19, 600x240')
    plt.xlim((i, i+60))
    plt.ylim((-10,10))
    plt.savefig('04-18-sky_'+str(i)+'-'+str(i+60)+'galex.png')
    plt.clf()


##########################################################
#g vs g-r
##########################################################
s1 = v2[np.where((v2['gl_sex'] > 0) & (v2['gl_sex'] < 40))]
s2 = v2[np.where((v2['gl_sex'] > 200) & (v2['gl_sex'] < 250))]
s3 = v2[np.where((v2['gl_sex'] > 250) & (v2['gl_sex'] < 300))]
s4 = v2[np.where((v2['gl_sex'] > 300) & (v2['gl_sex'] < 360))]

m = 8./1.5
b = 22 - m*0.5
#cut1 = np.where(((s1['g_AB']-s1['r_AB'])*m+b < s1['g_AB']) & (s1['g_AB'] > 18.5))
#cut2 = np.where(((s2['g_AB']-s2['r_AB'])*m+b < s2['g_AB']) & (s2['g_AB'] > 18.5))
#cut3 = np.where(((s3['g_AB']-s3['r_AB'])*m+b < s3['g_AB']) & (s3['g_AB'] > 18.5))
#cut4 = np.where(((s4['g_AB']-s4['r_AB'])*m+b < s4['g_AB']) & (s4['g_AB'] > 18.5))
vwd1 = vwd[np.where((vwd['gl_sex'] > 0) & (vwd['gl_sex'] < 40))]
vwd2 = vwd[np.where((vwd['gl_sex'] > 200) & (vwd['gl_sex'] < 250))]
vwd3 = vwd[np.where((vwd['gl_sex'] > 250) & (vwd['gl_sex'] < 300))]
vwd4 = vwd[np.where((vwd['gl_sex'] > 300) & (vwd['gl_sex'] < 360))]


fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
scatter_contour(s1['g_AB']-s1['r_AB'], s1['g_AB'], threshold=950, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0, 0])
axes[0, 0].scatter(vwd1['g_AB']-vwd1['r_AB'], vwd1['g_AB'], edgecolor='none', facecolor='blue')
axes[0, 0].axvline(x=0.59, c='black')
axes[0, 0].axhline(y=14.57, c='black')
axes[0, 0].set_ylabel(('g (ABmag)'))
axes[0, 0].set_xlim((-1, 3))
axes[0, 0].set_ylim((23, 12))
axes[0, 0].annotate('0<gl<40', xy=(1.75,14))

scatter_contour(s2['g_AB']-s2['r_AB'],s2['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0,1])
axes[0, 1].scatter(vwd2['g_AB']-vwd2['r_AB'],vwd2['g_AB'],edgecolor='none',facecolor='blue')
axes[0, 1].axvline(x=0.35, c='black')
axes[0, 1].axhline(y=15.1, c='black')
axes[0, 1].set_xlim((-1, 3))
axes[0, 1].set_ylim((23, 12))
axes[0,1].annotate('200<gl<250', xy=(1.75,14))

scatter_contour(s3['g_AB']-s3['r_AB'],s3['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 0])
axes[1, 0].scatter(vwd3['g_AB']-vwd3['r_AB'],vwd3['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 0].axvline(x=0.37, c='black')
axes[1, 0].axhline(y=14.5, c='black')
axes[1, 0].set_xlabel(('g - r (ABmag)'))
axes[1, 0].set_ylabel(('g (ABmag)'))
axes[1, 0].set_xlim((-1, 3))
axes[1, 0].set_ylim((23, 12))
axes[1, 0].annotate('250<gl<300', xy=(1.75,14))

scatter_contour(s4['g_AB']-s4['r_AB'],s4['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 1])
axes[1, 1].scatter(vwd4['g_AB']-vwd4['r_AB'],vwd4['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 1].axvline(x=0.52, c='black')
axes[1, 1].axhline(y=14.4, c='black')
axes[1, 1].set_xlabel(('g - r (ABmag)'))
axes[1, 1].set_xlim((-1, 3))
axes[1, 1].set_ylim((23, 12))
axes[1, 1].annotate('300<gl<360', xy=(1.75,14))
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle('vPHAS + SExtractor with WDCs')
plt.savefig('04-21-gvsgr_sex_vphas_allsky.png')


# Now with u-g < 1 cut
s1 = vcat[np.where((vcat['gl_sex'] > 0) & (vcat['gl_sex'] < 40) & (vcat['u_AB']-vcat['g_AB'] < 1))]
s2 = vcat[np.where((vcat['gl_sex'] > 200) & (vcat['gl_sex'] < 250) & (vcat['u_AB']-vcat['g_AB'] < 1))]
s3 = vcat[np.where((vcat['gl_sex'] > 250) & (vcat['gl_sex'] < 300) & (vcat['u_AB']-vcat['g_AB'] < 1))]
s4 = vcat[np.where((vcat['gl_sex'] > 300) & (vcat['gl_sex'] < 360) & (vcat['u_AB']-vcat['g_AB'] < 1))]

m = 8./1.5
b = 22 - m*0.5
cut1 = np.where(((s1['g_AB']-s1['r_AB'])*m+b < s1['g_AB']) & (s1['g_AB'] > 19.))
cut2 = np.where(((s2['g_AB']-s2['r_AB'])*m+b < s2['g_AB']) & (s2['g_AB'] > 19.))
cut3 = np.where(((s3['g_AB']-s3['r_AB'])*m+b < s3['g_AB']) & (s3['g_AB'] > 19.))
cut4 = np.where(((s4['g_AB']-s4['r_AB'])*m+b < s4['g_AB']) & (s4['g_AB'] > 19.))
wd1 = s1[cut1]
wd2 = s2[cut2]
wd3 = s3[cut3]
wd4 = s4[cut4]

scatter_contour(s1['g_AB']-s1['r_AB'],s1['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.scatter(wd1['g_AB']-wd1['r_AB'],wd1['g_AB'],edgecolor='none',facecolor='red')
plt.xlabel(('g - r (ABmag)'))
plt.ylabel(('g (ABmag)'))
plt.xlim((-1, 3))
plt.ylim((23, 12))
plt.title('VPHAS+SEx, 0 < gl < 40, u-g < 1, N$_{obj}$ = '+str(len(s1)))
plt.savefig('04-18-gvsgr_sex_vphas_0-40_ugcut.png')
plt.clf()

scatter_contour(s2['g_AB']-s2['r_AB'],s2['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.scatter(wd2['g_AB']-wd2['r_AB'],wd2['g_AB'],edgecolor='none',facecolor='red')
plt.xlabel(('g - r (ABmag)'))
plt.ylabel(('g (ABmag)'))
plt.xlim((-1, 3))
plt.ylim((23, 12))
plt.title('VPHAS+SEx, 200 < gl < 250, u-g < 1, N$_{obj}$ = '+str(len(s2)))
plt.savefig('04-18-gvsgr_sex_vphas_200-250_ugcut.png')
plt.clf()

scatter_contour(s3['g_AB']-s3['r_AB'],s3['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.scatter(wd3['g_AB']-wd3['r_AB'],wd3['g_AB'],edgecolor='none',facecolor='red')
plt.xlabel(('g - r (ABmag)'))
plt.ylabel(('g (ABmag)'))
plt.xlim((-1, 3))
plt.ylim((23, 12))
plt.title('VPHAS+SEx, 250 < gl < 300, u-g < 1, N$_{obj}$ = '+str(len(s3)))
plt.savefig('04-18-gvsgr_sex_vphas_250-300_ugcut.png')
plt.clf()

scatter_contour(s4['g_AB']-s4['r_AB'],s4['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.scatter(wd4['g_AB']-wd4['r_AB'],wd4['g_AB'],edgecolor='none',facecolor='red')
plt.xlabel(('g - r (ABmag)'))
plt.ylabel(('g (ABmag)'))
plt.xlim((-1, 3))
plt.ylim((23, 12))
plt.title('VPHAS+SEx, 300 < gl < 360, u-g < 1, N$_{obj}$ = '+str(len(s4)))
plt.savefig('04-18-gvsgr_sex_vphas_300-360_ugcut.png')
plt.clf()


###########################################################
# g-r vs nuv-g
###########################################################
s1 = v2[np.where((v2['gl_sex'] > 0) & (v2['gl_sex'] < 40))]
s2 = v2[np.where((v2['gl_sex'] > 200) & (v2['gl_sex'] < 250))]
s3 = v2[np.where((v2['gl_sex'] > 250) & (v2['gl_sex'] < 300))]
s4 = v2[np.where((v2['gl_sex'] > 300) & (v2['gl_sex'] < 360))]

pickles = Table.read('../picklemags_laphare_final.txt', format='ascii')
scatter_contour(s1['nuv_sex']-s1['g_AB'],s1['g_AB']-s1['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, NUV < 19, 0 < gl < 40, N$_{obj}$ = '+str(len(s1)))
plt.savefig('04-18-grvsnuvg_sex_vphas_0-40_nuvcut.png')
plt.clf()

scatter_contour(s2['nuv_sex']-s2['g_AB'],s2['g_AB']-s2['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, NUV < 19, 200 < gl < 250, N$_{obj}$ = '+str(len(s2)))
plt.savefig('04-18-grvsnuvg_sex_vphas_200-250_nuvcut.png')
plt.clf()

scatter_contour(s3['nuv_sex']-s3['g_AB'],s3['g_AB']-s3['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, NUV < 19, 250 < gl < 300, N$_{obj}$ = '+str(len(s3)))
plt.savefig('04-18-grvsnuvg_sex_vphas_250-300_nuvcut.png')
plt.clf()

scatter_contour(s4['nuv_sex']-s4['g_AB'],s4['g_AB']-s4['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, NUV < 19, 300 < gl < 360, N$_{obj}$ = '+str(len(s4)))
plt.savefig('04-18-grvsnuvg_sex_vphas_300-360_nuvcut.png')
plt.clf()

s1 = vcat[np.where((vcat['gl_sex'] > 0) & (vcat['gl_sex'] < 40))]
s2 = vcat[np.where((vcat['gl_sex'] > 200) & (vcat['gl_sex'] < 250))]
s3 = vcat[np.where((vcat['gl_sex'] > 250) & (vcat['gl_sex'] < 300))]
s4 = vcat[np.where((vcat['gl_sex'] > 300) & (vcat['gl_sex'] < 360))]

pickles = Table.read('../picklemags_laphare_final.txt', format='ascii')
scatter_contour(s1['nuv_sex']-s1['g_AB'],s1['g_AB']-s1['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, 0 < gl < 40, N$_{obj}$ = '+str(len(s1)))
plt.savefig('04-18-grvsnuvg_sex_vphas_0-40.png')
plt.clf()

scatter_contour(s2['nuv_sex']-s2['g_AB'],s2['g_AB']-s2['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, 200 < gl < 250, N$_{obj}$ = '+str(len(s2)))
plt.savefig('04-18-grvsnuvg_sex_vphas_200-250.png')
plt.clf()

scatter_contour(s3['nuv_sex']-s3['g_AB'],s3['g_AB']-s3['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, 250 < gl < 300, N$_{obj}$ = '+str(len(s3)))
plt.savefig('04-18-grvsnuvg_sex_vphas_250-300.png')
plt.clf()

scatter_contour(s4['nuv_sex']-s4['g_AB'],s4['g_AB']-s4['r_AB'],threshold=1050,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1))
plt.xlabel(('NUV - g(ABmag)'))
plt.ylabel(('g - r (ABmag)'))
plt.xlim((-4, 8))
plt.ylim((-1, 3))
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen')
plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.title('VPHAS+SEx, 300 < gl < 360, N$_{obj}$ = '+str(len(s4)))
plt.savefig('04-18-grvsnuvg_sex_vphas_300-360.png')
plt.clf()


a0 = gcat[np.where((gcat['gl_galex'] > 150) & (gcat['gl_galex'] < 155))]
a1 = gcat[np.where((gcat['gl_galex'] > 350) & (gcat['gl_galex'] < 355))]

m0, med0, std0 = [], [], []
m1, med1, std1 = [], [], []
mean0, median0, stdev0 = [], [], []
mean1, median1, stdev1 = [], [], []
mean2, median2, stdev2 = [], [], []
mean3, median3, stdev3 = [], [], []

for magrange in np.arange(11.5, 22, 0.5):
    mag0 = np.where((a0['nuv_galex'] > magrange) & (a0['nuv_galex'] < magrange+1))
    m0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
    std0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
    med0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

    mag1 = np.where((a1['nuv_galex'] > magrange) & (a1['nuv_galex'] < magrange+1))
    m1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
    std1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
    med1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

    mean0.append(np.mean((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
    median0.append(np.median((a0['nuv_sex']-a0['nuv_galex'])[mag0]))
    stdev0.append(np.std((a0['nuv_sex']-a0['nuv_galex'])[mag0]))

    mean1.append(np.mean((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
    median1.append(np.median((a1['nuv_sex']-a1['nuv_galex'])[mag1]))
    stdev1.append(np.std((a1['nuv_sex']-a1['nuv_galex'])[mag1]))

fig, (ax0, ax1) = plt.subplots(2, sharex=True, sharey=True)

scatter_contour(a0['nuv_galex'],a0['nuv_sex']-a0['nuv_galex'],threshold=500,log_counts=True,histogram2d_args=dict(bins=(30)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax0)
ax0.errorbar(np.arange(11.5, 22, 0.5), m0, yerr=std0, color='red', linewidth='2')
ax0.axhline(y=0, c='green')

scatter_contour(a1['nuv_galex'],a1['nuv_sex']-a1['nuv_galex'],threshold=500,log_counts=True,histogram2d_args=dict(bins=(30)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
ax1.errorbar(np.arange(11.5, 22, 0.5), m1, yerr=std1, color='red', linewidth='2')
ax1.axhline(y=0, c='green')
ax1.set_xlabel('NUV$_{GAIS}$')
ax0.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
ax1.set_ylabel('NUV$_{SEx}$ - NUV$_{GAIS}$')
ax0.set_xlim((12, 23))
ax0.set_ylim((-3, 2))
ax1.set_xlim((12, 23))
ax1.set_ylim((-3, 2))
ax0.annotate('gl = 150-155', xy=(13, -2))
ax1.annotate('gl = 350-355', xy=(13, -2))
fig.subplots_adjust(hspace=0)
plt.setp([lab.get_xticklabels() for lab in fig.axes[:-1]], visible=False)
plt.show()



##################################################################
# Find duplicates
##################################################################
skyrange = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']


a = Table.read('starcat_44.6mapweight.txt', format='ascii')
b = Table.read('starcat_45.5mapweight.txt', format='ascii')
agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
bgal = SkyCoord(b['gl']*u.deg, b['gb']*u.deg, frame='galactic')
aind, bind, angsep, ang3d = search_around_sky(agal, bgal, 1*u.arcsec)
if len(bind) != 0:
    b.remove_rows(bind)
alldata = vstack([a,b])

for region in range(0,len(skyrange)):
    a = Table.read('starcat_'+skyrange[region]+'mapweight.txt', format='ascii')
    try:
        b = Table.read('starcat_'+skyrange[region+1]+'mapweight.txt', format='ascii')
    except IndexError:
        ascii.write(alldata, 'starcat_all_nomatches_mapweight90-180.txt',format='basic')
    agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
    bgal = SkyCoord(b['gl']*u.deg, b['gb']*u.deg, frame='galactic')
    aind, bind, angsep, ang3d = search_around_sky(agal, bgal, 1*u.arcsec)
    if len(aind) != 0:
        a.remove_rows(aind)
    alldata = vstack([alldata, a])
    print skyrange[region]



a = Table.read('starcat_5mapweight.txt', format='ascii')
b = Table.read('starcat_1.4mapweight.txt', format='ascii')
agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
bgal = SkyCoord(b['gl']*u.deg, b['gb']*u.deg, frame='galactic')
aind, bind, angsep, ang3d = search_around_sky(agal, bgal, 2*u.arcsec)
a2 = a[aind]
b2 = b[bind]
plt.scatter(a['gl'], a['gb'], edgecolor='none')
plt.scatter(b['gl'], b['gb'], c='red', edgecolor='none')
plt.scatter(a2['gl'], a2['gb'], color='green', alpha=0.5)
plt.scatter(b2['gl'], b2['gb'], color='orange', alpha=0.5)



skyrange = ['1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']
alldata = Table.read('starcat_5mapweight.txt', format='ascii')
alldata = alldata[np.where((alldata['gl'] < (np.max(alldata['gl']) - 0.24)) & (alldata['nuv'] < 20))]

for region in skyrange:
    a = Table.read('starcat_'+region+'mapweight.txt', format='ascii')
    a = a[np.where((a['gl'] < (np.max(a['gl']) - 0.24)) & (a['nuv'] < 20))]
    alldata = vstack([alldata, a])

#################################################
#normalize histogram:
#################################################

area = 2000 #6118 # 360*20 - blank area
x, bins, p = plt.hist(vwd['gb_sex'], bins=50)
for item in p:
    item.set_height(item.get_height()/area)
plt.show()


# Histogram plot for sky map
h, xed, yed = np.histogram2d(a1['gl'],a1['gb'], bins=(np.linspace(0, 90, 1500), np.linspace(-10,10,1500)))
plt.imshow(h.T,vmin=0,vmax=20,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray, extent=[0,90,-10,10])
plt.show()

#################################################
# Subplots for each scan
#################################################
skyrange = ['5', '1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '9.5', '10.4', '11.3', '12.2', '14.0', '14.9', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '79.7', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '90.5', '91.4', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '110.3', '102.2', '103.1', '104.0', '104.9', '105.8', '106.7', '107.6', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3', '121.1', '122.9', '124.7', '125.6', '126.5', '127.4', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '223.7', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8', '241.7', '242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '273.2', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '283.1', '284.0', '285.8', '286.7', '288.5', '289.4', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '301.1', '302.0', '302.9', '303.8', '304.7', '305.6', '306.5', '308.3', '309.2', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '324.5', '325.4', '326.3', '327.2', '328.1', '329.0', '329.9', '331.7', '332.6', '333.5', '334.4', '335.3', '338.0', '338.9', '339.8', '341.6', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8', '358.7', '359.6']


tycho = Table(fits.open('../../tycho2.fits')[1].data)
tycho.remove_columns(('RAJ2000', 'DEJ2000', 'TYC1', 'TYC2', 'TYC3', 'pmRA', 'pmDE', 'BTmag', 'e_BTmag', 'VTmag', 'e_VTmag', 'HIP', 'RA_ICRS_', 'DE_ICRS_'))
galex = Table(fits.open('../../GALEXAIS.fits')[1].data)
galex.remove_columns(('fuv_mag', 'ra', 'dec'))
galexcut = np.where(galex['nuv_mag'] == -999.)
galex.remove_rows(galexcut)
tychogal = SkyCoord(tycho['Glon']*u.deg, tycho['Glat']*u.deg, frame='galactic')
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')

for region in skyrange:
    print region
    t1 = Table.read('starcat_'+region+'mapweight_fwhm.txt',format='ascii')
    t1.remove_columns(('X_IMAGE', 'Y_IMAGE', 'FLUX_AUTO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'x_new', 'y_new', 'FLUXERR_AUTO', 'FLUX_APER','ra', 'dec'))
    t1gal = SkyCoord(t1['gl']*u.deg, t1['gb']*u.deg, frame='galactic')
    t1tyind, tychoind, angsepty, ang3d = search_around_sky(t1gal, tychogal, 3*u.arcsec)
    t1ty = hstack([t1[t1tyind], Table(tycho[tychoind])])
    t1ty['angsep'] = angsepty
    t1ty.rename_column('gl', 'gl_sex')
    t1ty.rename_column('gb', 'gb_sex')
    t1ty.rename_column('Glon', 'gl_tycho')
    t1ty.rename_column('Glat', 'gb_tycho')
    
    if region == '5':
        for i in range(len(t1ty)):
            if t1ty['gl_tycho'][i] > 350:
                t1ty['gl_tycho'][i] = t1ty['gl_tycho'][i] - 360
    
    t1gaisind, galexind, angsepga, ang3d = search_around_sky(t1gal, galexgal, 3*u.arcsec)
    t1gais = hstack([t1[t1gaisind], Table(galex[galexind])])
    t1gais['angsep'] = angsepga
    t1gais.rename_column('nuv', 'nuv_sex')
    t1gais.rename_column('nuv_mag', 'nuv_galex')
    
    fig, axes = plt.subplots(2,2)
    axes[0,0].scatter((t1ty['gl_sex']-t1ty['gl_tycho'])*3600,(t1ty['gb_sex']-t1ty['gb_tycho'])* 3600, edgecolor='none', alpha=0.2)
    axes[0,0].axhline(y=0, c='red')
    axes[0,0].axvline(x=0, c='red')
    axes[0,0].set_xlim((-3, 3))
    axes[0,0].set_ylim((-3, 3))
    axes[0,0].set_title('gl (SEx-T2)')
    axes[0,0].set_ylabel('gb (SEx-T2)')
    
    scatter_contour(t1gais['nuv_galex'],t1gais['nuv_sex']-t1gais['nuv_galex'],threshold=500,    log_counts=True,histogram2d_args=dict(bins=(30)),plot_args=dict(color='k',markersize=1),    contour_args=dict(cmap=cm.gray), ax=axes[0,1])
    axes[0,1].axhline(y=0, c='red')
    axes[0,1].set_xlim((12, 25))
    axes[0,1].set_ylim((-2, 2))
    axes[0,1].set_ylabel('NUV (SEx-GAIS)')
    
    axes[1,0].hist(t1['nuv'], bins=50)
    axes[1,0].set_xlim((12,25))
    axes[1,0].set_xlabel('NUV (all SExtractor)')
    
    scatter_contour(t1['nuv'],t1['FWHM_IMAGE'],threshold=700,log_counts=True,histogram2d_args=dict  (bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes    [1,1])
    axes[1,1].axhline(y=0, c='red')
    axes[1,1].set_xlim((12,25))
    axes[1,1].set_ylim((-3, 12))
    axes[1,1].set_xlabel('NUV')
    axes[1,1].set_ylabel('FWHM')
    
    plt.suptitle('region = '+region)
    plt.savefig('region'+region+'matchplots.png')



# FWHM vs NUV
for region in range(0, 361, 5):
    a1 = a[np.where((a['gl'] > region) & (a['gl'] < region+5) & (a['FWHM_IMAGE'] > -10) & (a['FWHM_IMAGE'] < 100))]

    scatter_contour(a1['nuv'],a1['FWHM_IMAGE'],threshold=1500,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
    plt.xlim((14,25))
    plt.ylim((-3, 12))
    plt.axhline(y=0, c='red')
    plt.xlabel('NUV')
    plt.ylabel('FWHM')
    plt.title('gl = '+str(region)+'-'+str(region+5))
    plt.savefig('fwhmvsnuv'+str(region)+'-'+str(region+5)+'.png')
    plt.clf()