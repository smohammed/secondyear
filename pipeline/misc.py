######################################################################
# Reddening values 
######################################################################
NUV = 2.9720
B = 1.3172
V = 1.0
J = 0.2876
H = 0.1783
K = 0.1170
u = 1.5812
g = 1.1936
r = 0.8694
i = 0.6533

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
#sv = fits.open('sex_vphas.fits')[1].data
#vwd = Table.read('wds_sex_vphas.txt', format='ipac')
sv = fits.open('starcat_ps1_g10-20_09-05-17.fits')[1].data
pickles = Table.read('picklemags_laphare_final.txt', format='ascii')
#cmd = Table.read('cmdfiles/cmd_merged.txt', format='ascii')

#scatter_contour(sv['nuv']-sv['g_ps'],sv['g_ps']-sv['r_ps'],threshold=1100,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
plt.scatter(ps['nuv']-ps['g_ps'], ps['g_ps']-ps['r_ps'], color='k', s=1, alpha=0.06)

#plt.scatter(vwd['nuv_mag']-vwd['g_AB'],vwd['g_AB']-vwd['r_AB'], edgecolor='none', facecolor='blue', label='WDCs')
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='red', label='SED model', s=30)
#plt.plot(cmd['NUV']-cmd['g'],cmd['g']-cmd['r'],color='green', label='CMD')
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,0.6),2.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,1.1),6.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((2,-0.1),6,0.6,facecolor='lightblue',alpha=0.5,angle=5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-2,-0.75),3,0.65,facecolor='gray',alpha=0.5))
plt.arrow(3, -0.75, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.xlim((-3,9))
plt.ylim((-1,2))

plt.annotate('O', xy=(0, -0.9), size=20)
plt.annotate('B', xy=(1, -0.81), size=20)
plt.annotate('A', xy=(2.5, -0.5), size=20)
plt.annotate('F', xy=(3.7, -0.4), size=20)
plt.annotate('G', xy=(5, -0.2), size=20)
plt.annotate('K', xy=(6.2, -0.1), size=20)

plt.xlabel('NUV - g (ABmag)')
plt.ylabel('g - r (ABmag)')
plt.legend(scatterpoints=1, loc=4)
plt.title('vPHAS + SExtractor')
plt.show()

####################################################################
# NUV - g vs g - i
####################################################################
scatter_contour(v2['g_AB']-v2['i_AB'],v2['nuv_sex']-v2['g_AB'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
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
wdc1 = wdc[np.where(wdc['logg'] == 7)]
pickles = Table.read('../../picklemags_laphare_final.txt', format='ascii')

scatter_contour(cat['u_AB']-cat['g_AB'], cat['g_AB']-cat['r_AB'],threshold=1100,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
#plt.scatter(wd['u_AB']-wd['g_AB'], wd['g_AB']-wd['r_AB'],c='blue', edgecolor='none', s=30, label='WDCs')
#plt.scatter(var['u_AB']-var['g_AB'], var['g_AB']-var['r_AB'],c='purple', edgecolor='none', s=30, label='var')
#plt.plot(wdc1['u']-wdc1['g'], wdc1['g']-wdc1['r'],c='orange', label='DA WD')

plt.scatter(comb['u_AB']-comb['g_AB'], comb['g_AB']-comb['r_AB'], edgecolor='none', c='blue',s=20, label='vphas+gaia')
plt.scatter(pickles['u']-pickles['g'], pickles['g']-pickles['r'],c='red', edgecolor='none', s=30, label='SED model')

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
#plt.title('WDCs')
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


##########################################################
#g vs g-r
##########################################################
s1 = sv[np.where((sv['gl_sex'] > 0) & (sv['gl_sex'] < 40))]
s2 = sv[np.where((sv['gl_sex'] > 200) & (sv['gl_sex'] < 250))]
s3 = sv[np.where((sv['gl_sex'] > 250) & (sv['gl_sex'] < 300))]
s4 = sv[np.where((sv['gl_sex'] > 300) & (sv['gl_sex'] < 360))]
cmd = Table.read('cmdfiles/cmd_merged_zt.txt', format='ascii')

m = 8./1.5
b = 22 - m*0.5
cut1 = np.where(((s1['g_AB']-s1['r_AB'])*m+b < s1['g_AB']))
cut2 = np.where(((s2['g_AB']-s2['r_AB'])*m+b < s2['g_AB']))
cut3 = np.where(((s3['g_AB']-s3['r_AB'])*m+b < s3['g_AB']))
cut4 = np.where(((s4['g_AB']-s4['r_AB'])*m+b < s4['g_AB']))
vwd1 = s1[cut1]
vwd2 = s2[cut2]
vwd3 = s3[cut3]
vwd4 = s4[cut4]
#vwd1 = vwd[np.where((vwd['gl_sex'] > 0) & (vwd['gl_sex'] < 40))]
#vwd2 = vwd[np.where((vwd['gl_sex'] > 200) & (vwd['gl_sex'] < 250))]
#vwd3 = vwd[np.where((vwd['gl_sex'] > 250) & (vwd['gl_sex'] < 300))]
#vwd4 = vwd[np.where((vwd['gl_sex'] > 300) & (vwd['gl_sex'] < 360))]


fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
scatter_contour(s1['g_AB']-s1['r_AB'], s1['g_AB'], threshold=950, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0, 0])
axes[0, 0].scatter(vwd1['g_AB']-vwd1['r_AB'], vwd1['g_AB'], edgecolor='none', facecolor='blue')
axes[0, 0].plot(cmd['g']-cmd['r'], cmd['g'], c='green')
axes[0, 0].axvline(x=0.576, c='red')
axes[0, 0].axhline(y=14.272, c='red')
axes[0, 0].set_ylabel(('g (ABmag)'))
axes[0, 0].set_xlim((-1, 3))
axes[0, 0].set_ylim((23, 12))
axes[0, 0].annotate('0<gl<40', xy=(1.75,14))

scatter_contour(s2['g_AB']-s2['r_AB'],s2['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0,1])
axes[0, 1].plot(cmd['g']-cmd['r'], cmd['g'], c='green')

axes[0, 1].scatter(vwd2['g_AB']-vwd2['r_AB'],vwd2['g_AB'],edgecolor='none',facecolor='blue')
axes[0, 1].axvline(x=0.271, c='red')
axes[0, 1].axhline(y=14.69, c='red')
axes[0, 1].set_xlim((-1, 3))
axes[0, 1].set_ylim((23, 12))
axes[0,1].annotate('200<gl<250', xy=(1.75,14))

scatter_contour(s3['g_AB']-s3['r_AB'],s3['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 0])
axes[1, 0].scatter(vwd3['g_AB']-vwd3['r_AB'],vwd3['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 0].plot(cmd['g']-cmd['r'], cmd['g'], c='green')
axes[1, 0].axvline(x=0.369, c='red')
axes[1, 0].axhline(y=14.25, c='red')
axes[1, 0].set_xlabel(('g - r (ABmag)'))
axes[1, 0].set_ylabel(('g (ABmag)'))
axes[1, 0].set_xlim((-1, 3))
axes[1, 0].set_ylim((23, 12))
axes[1, 0].annotate('250<gl<300', xy=(1.75,14))

scatter_contour(s4['g_AB']-s4['r_AB'],s4['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 1])
axes[1, 1].scatter(vwd4['g_AB']-vwd4['r_AB'],vwd4['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 1].plot(cmd['g']-cmd['r'], cmd['g'], c='green')
axes[1, 1].axvline(x=0.504, c='red')
axes[1, 1].axhline(y=14.1, c='red')
axes[1, 1].set_xlabel(('g - r (ABmag)'))
axes[1, 1].set_xlim((-1, 3))
axes[1, 1].set_ylim((23, 12))
axes[1, 1].annotate('300<gl<360', xy=(1.75,14))
fig.subplots_adjust(wspace=0, hspace=0)
#fig.suptitle('11-22, SEx+vphas WDC cut')
plt.show()
#plt.savefig('11-22-gvsgr_sex_vphas.png')


# Now with u-g < 1 cut
s1 = sv[np.where((sv['gl_sex'] > 0) & (sv['gl_sex'] < 40) & (sv['u_AB']-sv['g_AB'] < 1))]
s2 = sv[np.where((sv['gl_sex'] > 200) & (sv['gl_sex'] < 250) & (sv['u_AB']-sv['g_AB'] < 1))]
s3 = sv[np.where((sv['gl_sex'] > 250) & (sv['gl_sex'] < 300) & (sv['u_AB']-sv['g_AB'] < 1))]
s4 = sv[np.where((sv['gl_sex'] > 300) & (sv['gl_sex'] < 360) & (sv['u_AB']-sv['g_AB'] < 1))]

m = 8./1.5
b = 22 - m*0.5
cut1 = np.where(((s1['g_AB']-s1['r_AB'])*m+b < s1['g_AB']) & (s1['g_AB'] > 19.))
cut2 = np.where(((s2['g_AB']-s2['r_AB'])*m+b < s2['g_AB']) & (s2['g_AB'] > 19.))
cut3 = np.where(((s3['g_AB']-s3['r_AB'])*m+b < s3['g_AB']) & (s3['g_AB'] > 19.))
cut4 = np.where(((s4['g_AB']-s4['r_AB'])*m+b < s4['g_AB']) & (s4['g_AB'] > 19.))
vwd1 = s1[cut1]
vwd2 = s2[cut2]
vwd3 = s3[cut3]
vwd4 = s4[cut4]


fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
scatter_contour(s1['g_AB']-s1['r_AB'], s1['g_AB'], threshold=950, log_counts=True, histogram2d_args=dict(bins=(40)), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0, 0])
axes[0, 0].scatter(vwd1['g_AB']-vwd1['r_AB'], vwd1['g_AB'], edgecolor='none', facecolor='blue')
axes[0, 0].axvline(x=0.576, c='red')
axes[0, 0].axhline(y=14.272, c='red')
axes[0, 0].set_ylabel(('g (ABmag)'))
axes[0, 0].set_xlim((-1, 3))
axes[0, 0].set_ylim((23, 12))
axes[0, 0].annotate('0<gl<40', xy=(1.75,14))

scatter_contour(s2['g_AB']-s2['r_AB'],s2['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[0,1])
axes[0, 1].scatter(vwd2['g_AB']-vwd2['r_AB'],vwd2['g_AB'],edgecolor='none',facecolor='blue')
axes[0, 1].axvline(x=0.271, c='red')
axes[0, 1].axhline(y=14.69, c='red')
axes[0, 1].set_xlim((-1, 3))
axes[0, 1].set_ylim((23, 12))
axes[0,1].annotate('200<gl<250', xy=(1.75,14))

scatter_contour(s3['g_AB']-s3['r_AB'],s3['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 0])
axes[1, 0].scatter(vwd3['g_AB']-vwd3['r_AB'],vwd3['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 0].axvline(x=0.369, c='red')
axes[1, 0].axhline(y=14.25, c='red')
axes[1, 0].set_xlabel(('g - r (ABmag)'))
axes[1, 0].set_ylabel(('g (ABmag)'))
axes[1, 0].set_xlim((-1, 3))
axes[1, 0].set_ylim((23, 12))
axes[1, 0].annotate('250<gl<300', xy=(1.75,14))

scatter_contour(s4['g_AB']-s4['r_AB'],s4['g_AB'],threshold=950,log_counts=True,histogram2d_args=dict(bins=(40)),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=axes[1, 1])
axes[1, 1].scatter(vwd4['g_AB']-vwd4['r_AB'],vwd4['g_AB'],edgecolor='none',facecolor='blue')
axes[1, 1].axvline(x=0.504, c='red')
axes[1, 1].axhline(y=14.1, c='red')
axes[1, 1].set_xlabel(('g - r (ABmag)'))
axes[1, 1].set_xlim((-1, 3))
axes[1, 1].set_ylim((23, 12))
axes[1, 1].annotate('300<gl<360', xy=(1.75,14))
fig.subplots_adjust(wspace=0, hspace=0)
fig.suptitle('11-22, SEx+vphas, WDC cut, u-g < 1')
#plt.show()
plt.savefig('11-22-gvsgr_sex_vphas_ugcut.png')


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


#################################################
# Combine all starcat files
#################################################
incscans = ['14.9', '79.7', '90.5', '91.4', '103.1', '104.0', '122.9', '127.4', '223.7', '273.2', '283.1', '289.4', '306.5', '309.2', '324.5', '329.9', '338.0', '339.8', '342.5', '343.4', '345.2', '348.8', '349.7', '350.6', '351.5', '352.4', '353.3', '354.2', '355.1', '356.0', '357.8']

alldata = Table.read('starcat_9.5mapweight_fwhm.txt', format='ascii')
alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5))]

for region in scans:
    a = Table.read('starcat_'+region+'mapweight_fec_fwhm.txt', format='ascii')
    #a = a[np.where((a['FWHM_IMAGE'] < 10) & (a['FWHM_IMAGE'] > 3.5) & (a['nuv'] < 21))]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, '../../../starcat_incscans_mapweight_fwhm_pscans_allnuv.txt', format='ipac')


sky1 = ['1.4', '2.3', '3.2', '4.1', '5.0', '5.9', '6.8', '8.6', '10.4', '11.3', '12.2', '14.0', '15.8', '16.7', '17.6', '18.5', '19.4', '20.3', '21.2', '22.1', '23.0', '23.9', '24.8', '25.7', '28.4', '29.3', '30.2', '31.1', '32.0', '32.9', '33.8', '34.7', '35.6', '39.2', '42.8', '43.7', '44.6', '44.6', '45.5', '46.4', '47.3', '48.2', '49.1', '50.0', '67.1', '68.9', '71.6', '74.3', '75.2', '76.1', '77.0', '77.9', '78.8', '80.6', '81.5', '82.4', '83.3', '87.8', '88.7', '89.6', '92.3', '93.2', '94.1', '95.0', '95.9', '96.8', '97.7', '98.6', '99.5', '100.4', '101.3', '102.2', '104.9', '105.8', '106.7', '107.6', '110.3', '111.2', '112.1', '113.0', '113.9', '114.8', '119.3']

alldata = Table.read('starcat_5mapweight_fwhm.txt', format='ascii')
alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5) & (alldata['nuv'] < 21))]

for region in sky1:
    a = Table.read('starcat_'+region+'mapweight_fwhm.txt', format='ascii')
    a = a[np.where((a['FWHM_IMAGE'] < 10) & (a['FWHM_IMAGE'] > 3.5) & (a['nuv'] < 21))]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, '../../../starcat_1-120_mapweight_fwhm_pscans_allnuv.txt', format='ipac')


sky2 = ['124.7', '125.6', '126.5', '128.3', '129.2', '130.1', '131.0', '131.9', '132.8', '133.7', '134.6', '135.5', '136.4', '137.3', '138.2', '139.1', '140.0', '140.9', '141.8', '143.6', '144.5', '148.1', '149.0', '149.9', '150.8', '151.7', '152.6', '153.5', '145.4', '155.3', '156.2', '157.1', '158.0', '160.7', '161.6', '163.4', '167.0', '167.9', '172.4', '173.3', '174.2', '175.1', '176.0', '176.9', '177.8', '178.7', '179.6', '180.5', '183.2', '185.0', '190.4', '191.3', '197.6', '198.5', '200.3', '201.2', '203.0', '203.9', '205.7', '206.6', '207.5', '208.4', '209.3', '210.2', '211.1', '212.0', '212.9', '213.8', '214.7', '215.6', '216.5', '217.4', '218.3', '219.2', '220.1', '221.0', '221.9', '222.8', '224.6', '225.5', '226.4', '228.2', '229.1', '230.0', '230.9', '231.8', '234.5', '235.4', '236.3', '237.2', '238.1', '239.0', '239.9', '240.8']

alldata = Table.read('starcat_121.1mapweight_fwhm.txt', format='ascii')
alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5) & (alldata['nuv'] < 21))]

for region in sky2:
    a = Table.read('starcat_'+region+'mapweight_fwhm.txt', format='ascii')
    a = a[np.where((a['FWHM_IMAGE'] < 10) & (a['FWHM_IMAGE'] > 3.5) & (a['nuv'] < 21))]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, '../../../starcat120-240_mapweight_fwhm_pscans_allnuv.txt', format='ipac')

sky3 =  ['242.6', '243.5', '244.4', '245.3', '246.2', '247.1', '248.0', '248.9', '249.8', '250.7', '251.6', '252.5', '253.4', '254.3', '255.2', '256.1', '257.0', '258.8', '259.7', '260.6', '261.5', '263.3', '264.2', '265.1', '266.0', '266.9', '268.7', '269.6', '270.5', '271.4', '272.3', '274.1', '275.0', '275.9', '276.8', '278.6', '279.5', '281.3', '284.0', '285.8', '286.7', '288.5', '290.3', '291.2', '292.1', '293.0', '293.9', '295.7', '297.5', '298.4', '302.0', '302.9', '303.8', '304.7', '305.6', '308.3', '310.1', '315.5', '316.4', '317.3', '318.2', '319.1', '320.0', '320.9', '321.8', '322.7', '323.6', '325.4', '326.3', '327.2', '328.1', '329.0', '331.7', '332.6', '333.5', '334.4', '335.3', '338.9', '341.6', '358.7', '359.6']

alldata = Table.read('starcat_241.7mapweight_fwhm.txt', format='ascii')
alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5) & (alldata['nuv'] < 21))]

for region in sky3:
    a = Table.read('starcat_'+region+'mapweight_fwhm.txt', format='ascii')
    a = a[np.where((a['FWHM_IMAGE'] < 10) & (a['FWHM_IMAGE'] > 3.5) & (a['nuv'] < 21))]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, '../../../starcat240-360_mapweight_fwhm_pscans_allnuv.txt', format='ipac')


#################################################
#normalize histogram:
#################################################
area = 100
x, bins, p = plt.hist(np.log10(sg['dist']), range=[0,4], bins=15, label='GAIS dist/100')
for item in p:
    item.set_height(item.get_height()/area)
plt.show()


# Histogram plot for sky map
h, xed, yed = np.histogram2d(a1['gl'],a1['gb'], bins=(np.linspace(0, 90, 1500), np.linspace(-10,10,1500)))
plt.imshow(h.T,vmin=0,vmax=20,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray, extent=[0,90,-10,10])
plt.show()


sky = ['329.0']
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

#################################################
# Combine all low/hi fec files
#################################################
files = np.loadtxt('lowfec.txt', dtype='str')
alldata = Table()
for region in files:
    print region
    a = Table.read(region, format='ascii')
    a = a[np.where((a['FWHM_IMAGE'] < 10) & (a['FWHM_IMAGE'] > 3.5))]
    alldata = vstack([alldata, a])

ascii.write(alldata, 'lowfec_all.txt', format='basic')


#################################################
# Match low/hi fec with catalogs
#################################################
galex = Table(fits.open('../../GALEXAIS.fits')[1].data)
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')

a = Table.read('starcat_'+reg+'mapweight_fwhm.txt',format='ascii')
agal = SkyCoord(a['gl']*u.deg, a['gb']*u.deg, frame='galactic')
aind, galexind, angsep, ang3d = search_around_sky(agal, galexgal, 3*u.arcsec)
a2 = a[aind]
plt.scatter(a2['gl'],a2['gb']), plt.show()


galex = Table(fits.open('../GALEXAIS.fits')[1].data)
galexgal = SkyCoord(galex['glon']*u.deg, galex['glat']*u.deg, frame='galactic')
high = fits.open('../highfec/sex_highfec.fits')[1].data
low = fits.open('../lowfec/lowfec_all.fits')[1].data
higal = SkyCoord(high['gl']*u.deg, high['gb']*u.deg, frame='galactic')
lowgal = SkyCoord(low['gl']*u.deg, low['gb']*u.deg, frame='galactic')
hiind, galhiind, anghi, ang3d = search_around_sky(higal, galexgal, 3*u.arcsec)
lowind, gallowind, anglow, ang3d = search_around_sky(lowgal, galexgal, 3*u.arcsec)

hi2 = Table(high[hiind])
ghi2 = Table(galex[galhiind])
hi2 = hstack([hi2, ghi2])

lo2 = Table(low[lowind])
glo2 = Table(galex[gallowind])
lo2 = hstack([lo2, glo2])

lo2.rename_column('nuv', 'nuv_sex')
lo2.rename_column('gl', 'gl_sex')
lo2.rename_column('gb', 'gb_sex')
lo2.rename_column('ra_1', 'ra_sex')
lo2.rename_column('dec_1', 'dec_sex')
lo2.rename_column('nuv_mag', 'nuv_galex')
lo2.rename_column('ra_2', 'ra_galex')
lo2.rename_column('dec_2', 'dec_galex')
lo2.rename_column('glon', 'gl_galex')
lo2.rename_column('glat', 'gb_galex')
lo2cut = np.where(lo2['nuv_galex'] == -999.)
lo2.remove_rows(lo2cut)

hi2.rename_column('nuv', 'nuv_sex')
hi2.rename_column('gl', 'gl_sex')
hi2.rename_column('gb', 'gb_sex')
hi2.rename_column('ra_1', 'ra_sex')
hi2.rename_column('dec_1', 'dec_sex')
hi2.rename_column('nuv_mag', 'nuv_galex')
hi2.rename_column('ra_2', 'ra_galex')
hi2.rename_column('dec_2', 'dec_galex')
hi2.rename_column('glon', 'gl_galex')
hi2.rename_column('glat', 'gb_galex')
hi2cut = np.where(hi2['nuv_galex'] == -999.)
hi2.remove_rows(hi2cut)

plt.scatter(lo2['nuv_galex'], lo2['nuv_sex']-lo2['nuv_galex'], edgecolor='none')
plt.scatter(hi2['nuv_galex'], hi2['nuv_sex']-hi2['nuv_galex'], edgecolor='none', c='red')


high = fits.open('../highfec/sex_highfec.fits')[1].data
low = fits.open('../lowfec/lowfec_all.fits')[1].data
higal = SkyCoord(high['gl']*u.deg, high['gb']*u.deg, frame='galactic')
lowgal = SkyCoord(low['gl']*u.deg, low['gb']*u.deg, frame='galactic')

vphas = fits.open('../vphas_allg.fits')[1].data
vpgal = SkyCoord(vphas['RAJ2000']*u.deg,vphas['DEJ2000']*u.deg).galactic


hiind, vphiind, anghi, ang3d = search_around_sky(higal, vpgal, 3*u.arcsec)
lowind, vplowind, anglow, ang3d = search_around_sky(lowgal, vpgal, 3*u.arcsec)

hi2 = Table(high[hiind])
ghi2 = Table(vphas[vphiind])
hi2 = hstack([hi2, ghi2])

lo2 = Table(low[lowind])
glo2 = Table(vphas[vplowind])
lo2 = hstack([lo2, glo2])

hi2['angsep'] = anghi
lo2['angsep'] = anglow

lo2.rename_column('nuv', 'nuv_sex')
lo2.rename_column('gl', 'gl_sex')
lo2.rename_column('gb', 'gb_sex')
lo2.rename_column('ra', 'ra_sex')
lo2.rename_column('dec', 'dec_sex')
lo2.rename_column('RAJ2000', 'ra_vphas')
lo2.rename_column('DEJ2000', 'dec_vphas')

hi2.rename_column('nuv', 'nuv_sex')
hi2.rename_column('gl', 'gl_sex')
hi2.rename_column('gb', 'gb_sex')
hi2.rename_column('ra', 'ra_sex')
hi2.rename_column('dec', 'dec_sex')
hi2.rename_column('RAJ2000', 'ra_vphas')
hi2.rename_column('DEJ2000', 'dec_vphas')

lo2 = lo2[np.unique(lo2['ra_vphas'], return_index=True)[1]]
hi2 = hi2[np.unique(hi2['ra_vphas'], return_index=True)[1]]

pickles = Table.read('../../picklemags_laphare_final.txt', format='ascii')


plt.arrow(6, -0.5, 2.972-1.1838, 1.1838-0.8664, head_length=0.05, head_width=0.02, color='red')
plt.scatter(pickles['nuv']-pickles['g'],pickles['g']-pickles['r'],color='darkgreen', label='pickles', s=30)
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,0.6),2.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-0.5,1.1),6.5,0.5,facecolor='yellow',alpha=0.5))
plt.gca().add_patch(matplotlib.patches.Rectangle((2,0.1),6,0.6,facecolor='lightblue',alpha=0.5,angle=5))
plt.gca().add_patch(matplotlib.patches.Rectangle((-2,-0.75),3,0.65,facecolor='gray',alpha=0.5))

comb.rename_column('nuv', 'nuv_sex')
comb.rename_column('gl', 'gl_sex')
comb.rename_column('gb', 'gb_sex')
comb.rename_column('ra_1', 'ra_sex')
comb.rename_column('dec_1', 'dec_sex')
comb.rename_column('nuv_mag', 'nuv_galex')
comb.rename_column('ra_2', 'ra_galex')
comb.rename_column('dec_2', 'dec_galex')
comb.rename_column('glon', 'gl_galex')
comb.rename_column('glat', 'gb_galex')
combcut = np.where(comb['nuv_galex'] == -999.)
comb.remove_rows(combcut)

#################################################
# Easy sky plot
#################################################
im = fits.open('starcat_all_mapweight_skymaponly.fits')[1].data
im1 = im[np.where((im['gl'] > 0) & (im['gl'] < 90))]
im2 = im[np.where((im['gl'] > 90) & (im['gl'] < 180))]
im3 = im[np.where((im['gl'] > 180) & (im['gl'] < 270))]
im4 = im[np.where((im['gl'] > 270) & (im['gl'] < 360))]
plt.hist2d(im1['gl'], im1['gb'], bins=[2000,400], cmap=cm.gray, vmin=0, vmax=15), plt.show()

#################################################
# WD cooling curve plots
#################################################
wdcda = Table.read('wdcool_DA.txt',format='ascii')
wdcdb = Table.read('wdcool_DB.txt',format='ascii')
obj = Table.read('Dunmaps/fwhm/wds_vstar_sex_mast_vphas_2mass.txt',format='ascii')
wd = obj[np.where(obj['type'] == 'wd')]


logg = np.arange(7, 9.1, 0.5)
for loggval in logg:
    wdcda1 = wdcda[np.where(wdcda['logg'] == loggval)]
    wdcdb1 = wdcdb[np.where(wdcdb['logg'] == loggval)]

    plt.scatter(wd['u_AB']-wd['g_AB'], wd['g_AB']-wd['r_AB'],c='blue', edgecolor='none', s=30, label='wd')
    plt.plot(wdcda1['u'] - wdcda1['g'], wdcda1['g']-wdcda1['r'], c='orange', linewidth='3', label='DA')
    plt.plot(wdcdb1['u'] - wdcdb1['g'], wdcdb1['g']-wdcdb1['r'], c='red', linewidth='3', label='DB')
    plt.xlabel('u - g (ABmag)')
    plt.ylabel('g - r (ABmag)')
    plt.xlim((-1, 1))
    plt.ylim((-1, 0.5))
    plt.legend(scatterpoints=1, loc=2)
    plt.title('WDCs, logg = '+str(loggval))
    plt.savefig('09-12-wd_dadb_curve_logg'+str(loggval)+'.png')
    plt.clf()
    #plt.show()


plt.scatter(wd['gl_sex'], wd['gb_sex'], edgecolor='none', label='WD')
plt.scatter(var['gl_sex'], var['gb_sex'], edgecolor='none', facecolor='red', label='var')
plt.xlabel('gl')
plt.ylabel('gb')
plt.title('WD + var dets in sextractor')


#################################################
# HR diagram with extinction colored by gl
#################################################
#sg = fits.open('../sex-gaia-dust.fits')[1].data
sg = fits.open('sex_gaia_dust_interp.fits')[1].data

sg1 = sg[np.where(sg['dist'] < 100)]
sg2 = sg[np.where((sg['dist'] > 100) & (sg['dist'] < 300))]
sg3 = sg[np.where((sg['dist'] > 300) & (sg['dist'] < 600))]
sg4 = sg[np.where((sg['dist'] > 600) & (sg['dist'] < 1000))]
sg5 = sg[np.where((sg['dist'] > 1000) & (sg['dist'] < 3000))]
sg6 = sg[np.where(sg['dist'] > 3000)]


# Extinction values from Schlafly & Finkbeiner 2011
scatter_contour(sg['nuv_mag']-sg['ebv']*7.76-sg['phot_g_mean_mag']-3.303, sg['Mg']-sg['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray))
plt.xlim((-6, 8))
plt.ylim((8, -6))
plt.xlabel('(NUV - E$_{B-V}$ * 7.76) - (g - E$_{B-V}$ * 3.303)')
plt.ylabel('Mg - E$_{B-V}$ * 3.303')
plt.title('S+G matches, Ext from S&F11')

#sg = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['dist'] < 3000))]
sg = sgcat[np.where((sgcat['dist'] > 3000))]

# By gl
sg1 = sg[np.where((sg['gl_sex'] > 0) & (sg['gl_sex'] < 90))]
sg2 = sg[np.where((sg['gl_sex'] > 90) & (sg['gl_sex'] < 180))]
sg3 = sg[np.where((sg['gl_sex'] > 180) & (sg['gl_sex'] < 270))]
sg4 = sg[np.where((sg['gl_sex'] > 270) & (sg['gl_sex'] < 360))]

fig = plt.figure()
plt.subplot(2, 2, 1)
plt.scatter(sg1['nuv_mag']-sg1['ebv']*7.76 - (sg1['phot_g_mean_mag']-3.303), sg1['Mg']-sg1['ebv']*3.303, edgecolor='none', c=sg1['ebv'], vmin=0, vmax=90, alpha=0.3)
plt.ylabel('Mg - E$_{B-V}$ * 3.303')
plt.xlim((-5, 17))
plt.ylim((10, -10))
plt.colorbar().set_label('gl (0 to 90)')

plt.subplot(2, 2, 2)
plt.scatter(sg2['nuv_mag']-sg2['ebv']*7.76 - (sg2['phot_g_mean_mag']-3.303), sg2['Mg']-sg2['ebv']*3.303, edgecolor='none', c=sg2['ebv'], vmin=90, vmax=180, alpha=0.3)
plt.xlim((-5, 17))
plt.ylim((10, -10))
plt.colorbar().set_label('gl (90 to 180)')

plt.subplot(2, 2, 3)
plt.scatter(sg3['nuv_mag']-sg3['ebv']*7.76 - (sg3['phot_g_mean_mag']-3.303), sg3['Mg']-sg3['ebv']*3.303, edgecolor='none', c=sg3['ebv'], vmin=180, vmax=270, alpha=0.3)
plt.xlabel('(NUV - E$_{B-V}$ * 7.76) - (g - E$_{B-V}$ * 3.303)')
plt.ylabel('Mg - E$_{B-V}$ * 3.303')
plt.xlim((-5, 17))
plt.ylim((10, -10))
plt.colorbar().set_label('gl (180 to 270)')

plt.subplot(2, 2, 4)
plt.scatter(sg4['nuv_mag']-sg4['ebv']*7.76 - (sg4['phot_g_mean_mag']-3.303), sg4['Mg']-sg4['ebv']*3.303, edgecolor='none', c=sg4['ebv'], vmin=270, vmax=360, alpha=0.3)
plt.xlabel('(NUV - E$_{B-V}$ * 7.76) - (g - E$_{B-V}$ * 3.303)')
plt.xlim((-5, 17))
plt.ylim((10, -10))
plt.colorbar().set_label('gl (270 to 360)')

plt.subplots_adjust(wspace=0.1, hspace=0)
fig.suptitle('S+G matches by gl, D > 3000 pc')

#######################################################
# HR diagram with extinction by dist and 45 deg slices
#######################################################
sgcat = fits.open('sex_gaia_dust_interp.fits')[1].data
#sgcat = fits.open('gais_gaia_dust.fits')[1].data
sga = sgcat[np.where((sgcat['dist'] > 0) & (sgcat['dist'] < 100) & (sgcat['ebv'] > 0))]
sgb = sgcat[np.where((sgcat['dist'] > 100) & (sgcat['dist'] < 300) & (sgcat['ebv'] > 0))]
sgc = sgcat[np.where((sgcat['dist'] > 300) & (sgcat['dist'] < 600) & (sgcat['ebv'] > 0))]
sgd = sgcat[np.where((sgcat['dist'] > 600) & (sgcat['dist'] < 1000) & (sgcat['ebv'] > 0))]
sge = sgcat[np.where((sgcat['dist'] > 1000) & (sgcat['dist'] < 3000) & (sgcat['ebv'] > 0))]
sgf = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['ebv'] > 0))]
cmd = Table.read('cmdfiles/cmd_merged.txt', format='ascii')

sg = [sga, sgb, sgc, sgd, sge]
distrange = ['0 < D < 100 pc', '100 < D < 300 pc', '300 < D < 600 pc', '600 < D < 1000 pc', 'D > 1000 pc']
distfilename = ['0-100', '100-300', '300-600', '600-1000', '1000']

horlinelim = [6.5, 6, 4.5, 3.5, 2.5]
verlinelim = [5, 4, 3.3, 3, 3]

for i in range(len(sg)):
    # By gl
    sg1 = sg[i][np.where((sg[i]['gl_sex'] > 0) & (sg[i]['gl_sex'] < 45))]
    sg2 = sg[i][np.where((sg[i]['gl_sex'] > 45) & (sg[i]['gl_sex'] < 90))]
    sg3 = sg[i][np.where((sg[i]['gl_sex'] > 90) & (sg[i]['gl_sex'] < 135))]
    sg4 = sg[i][np.where((sg[i]['gl_sex'] > 135) & (sg[i]['gl_sex'] < 180))]
    sg5 = sg[i][np.where((sg[i]['gl_sex'] > 180) & (sg[i]['gl_sex'] < 225))]
    sg6 = sg[i][np.where((sg[i]['gl_sex'] > 225) & (sg[i]['gl_sex'] < 270))]
    sg7 = sg[i][np.where((sg[i]['gl_sex'] > 270) & (sg[i]['gl_sex'] < 315))]
    sg8 = sg[i][np.where((sg[i]['gl_sex'] > 315) & (sg[i]['gl_sex'] < 360))]


    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True)
    scatter_contour((sg1['nuv_mag']-sg1['ebv']*7.76)-(sg1['phot_g_mean_mag']-sg1['ebv']*3.303), sg1['Mg']-sg1['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
    ax1.axhline(y=horlinelim[i], color='blue')
    ax1.axvline(x=verlinelim[i], color='blue')
    ax1.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg2['nuv_mag']-sg2['ebv']*7.76)-(sg2['phot_g_mean_mag']-sg2['ebv']*3.303), sg2['Mg']-sg2['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax2)
    ax2.axhline(y=horlinelim[i], color='blue')
    ax2.axvline(x=verlinelim[i], color='blue')
    ax2.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg3['nuv_mag']-sg3['ebv']*7.76)-(sg3['phot_g_mean_mag']-sg3['ebv']*3.303), sg3['Mg']-sg3['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax3)
    ax3.axhline(y=horlinelim[i], color='blue')
    ax3.axvline(x=verlinelim[i], color='blue')
    ax3.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg4['nuv_mag']-sg4['ebv']*7.76)-(sg4['phot_g_mean_mag']-sg4['ebv']*3.303), sg4['Mg']-sg4['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
    ax4.axhline(y=horlinelim[i], color='blue')
    ax4.axvline(x=verlinelim[i], color='blue')
    ax4.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg5['nuv_mag']-sg5['ebv']*7.76)-(sg5['phot_g_mean_mag']-sg5['ebv']*3.303), sg5['Mg']-sg5['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax5)
    ax5.axhline(y=horlinelim[i], color='blue')
    ax5.axvline(x=verlinelim[i], color='blue')
    ax5.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg6['nuv_mag']-sg6['ebv']*7.76)-(sg6['phot_g_mean_mag']-sg6['ebv']*3.303), sg6['Mg']-sg6['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax6)
    ax6.axhline(y=horlinelim[i], color='blue')
    ax6.axvline(x=verlinelim[i], color='blue')
    ax6.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg7['nuv_mag']-sg7['ebv']*7.76)-(sg7['phot_g_mean_mag']-sg7['ebv']*3.303), sg7['Mg']-sg7['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax7)
    ax7.axhline(y=horlinelim[i], color='blue')
    ax7.axvline(x=verlinelim[i], color='blue')
    ax7.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour((sg8['nuv_mag']-sg8['ebv']*7.76)-(sg8['phot_g_mean_mag']-sg8['ebv']*3.303), sg8['Mg']-sg8['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax8)
    ax8.axhline(y=horlinelim[i], color='blue')
    ax8.axvline(x=verlinelim[i], color='blue')
    ax8.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    ax1.set_xlim((-4, 11))
    ax1.set_ylim((8, -6))
    ax2.set_xlim((-4, 11))
    ax2.set_ylim((8, -6))
    ax3.set_xlim((-4, 11))
    ax3.set_ylim((8, -6))
    ax4.set_xlim((-4, 11))
    ax4.set_ylim((8, -6))
    ax5.set_xlim((-4, 11))
    ax5.set_ylim((8, -6))
    ax6.set_xlim((-4, 11))
    ax6.set_ylim((8, -6))
    ax7.set_xlim((-4, 11))
    ax7.set_ylim((8, -6))
    ax8.set_xlim((-4, 11))
    ax8.set_ylim((8, -6))
    ax1.annotate('gl = 0-45', xy=(-3.5, 7))
    ax2.annotate('gl = 45-90', xy=(-3.5, 7))
    ax3.annotate('gl = 90-135', xy=(-3.5, 7))
    ax4.annotate('gl = 135-180', xy=(-3.5, 7))
    ax5.annotate('gl = 180-225', xy=(-3.5, 7))
    ax6.annotate('gl = 225-270', xy=(-3.5, 7))
    ax7.annotate('gl = 270-315', xy=(-3.5, 7))
    ax8.annotate('gl = 315-360', xy=(-3.5, 7))
    fig.text(0.5, 0.04, '(NUV - E$_{B-V}$ * 7.76) - (G - E$_{B-V}$ * 3.303)', ha='center')
    fig.text(0.04, 0.5, 'MG - E$_{B-V}$ * 3.303', va='center', rotation='vertical')
    plt.suptitle('SEx+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', iso t=1e9yr, Z=0.0152')
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig('11-29-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'pc_ext_interp.png')
    plt.clf()
    #plt.show()


#######################################################
# HR diagram WITHOUT extinction by dist and 45 deg slices
#######################################################
sgcat = fits.open('sex_gaia_dust_interp.fits')[1].data
sga = sgcat[np.where((sgcat['dist'] > 0) & (sgcat['dist'] < 100) & (sgcat['ebv'] > 0))]
sgb = sgcat[np.where((sgcat['dist'] > 100) & (sgcat['dist'] < 300) & (sgcat['ebv'] > 0))]
sgc = sgcat[np.where((sgcat['dist'] > 300) & (sgcat['dist'] < 600) & (sgcat['ebv'] > 0))]
sgd = sgcat[np.where((sgcat['dist'] > 600) & (sgcat['dist'] < 1000) & (sgcat['ebv'] > 0))]
sge = sgcat[np.where((sgcat['dist'] > 1000) & (sgcat['dist'] < 3000) & (sgcat['ebv'] > 0))]
sgf = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['ebv'] > 0))]
cmd = Table.read('cmdfiles/cmd_merged.txt', format='ascii')

sg = [sga, sgb, sgc, sgd, sge]
distrange = ['0 < D < 100 pc', '100 < D < 300 pc', '300 < D < 600 pc', '600 < D < 1000 pc', 'D > 1000 pc']
distfilename = ['0-100', '100-300', '300-600', '600-1000', '1000']

horlinelim = [6.5, 6, 4.5, 3.5, 2.5]
verlinelim = [5, 4, 3.3, 3, 3]

for i in range(len(sg)):
    # By gl
    sg1 = sg[i][np.where((sg[i]['gl_sex'] > 0) & (sg[i]['gl_sex'] < 45))]
    sg2 = sg[i][np.where((sg[i]['gl_sex'] > 45) & (sg[i]['gl_sex'] < 90))]
    sg3 = sg[i][np.where((sg[i]['gl_sex'] > 90) & (sg[i]['gl_sex'] < 135))]
    sg4 = sg[i][np.where((sg[i]['gl_sex'] > 135) & (sg[i]['gl_sex'] < 180))]
    sg5 = sg[i][np.where((sg[i]['gl_sex'] > 180) & (sg[i]['gl_sex'] < 225))]
    sg6 = sg[i][np.where((sg[i]['gl_sex'] > 225) & (sg[i]['gl_sex'] < 270))]
    sg7 = sg[i][np.where((sg[i]['gl_sex'] > 270) & (sg[i]['gl_sex'] < 315))]
    sg8 = sg[i][np.where((sg[i]['gl_sex'] > 315) & (sg[i]['gl_sex'] < 360))]


    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True)
    scatter_contour(sg1['nuv_mag']-sg1['phot_g_mean_mag'], sg1['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
    ax1.axhline(y=horlinelim[i], color='blue')
    ax1.axvline(x=verlinelim[i], color='blue')
    ax1.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg2['nuv_mag']-sg2['phot_g_mean_mag'], sg2['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax2)
    ax2.axhline(y=horlinelim[i], color='blue')
    ax2.axvline(x=verlinelim[i], color='blue')
    ax2.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg3['nuv_mag']-sg3['phot_g_mean_mag'], sg3['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax3)
    ax3.axhline(y=horlinelim[i], color='blue')
    ax3.axvline(x=verlinelim[i], color='blue')
    ax3.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg4['nuv_mag']-sg4['phot_g_mean_mag'], sg4['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
    ax4.axhline(y=horlinelim[i], color='blue')
    ax4.axvline(x=verlinelim[i], color='blue')
    ax4.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg5['nuv_mag']-sg5['phot_g_mean_mag'], sg5['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax5)
    ax5.axhline(y=horlinelim[i], color='blue')
    ax5.axvline(x=verlinelim[i], color='blue')
    ax5.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg6['nuv_mag']-sg6['phot_g_mean_mag'], sg6['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax6)
    ax6.axhline(y=horlinelim[i], color='blue')
    ax6.axvline(x=verlinelim[i], color='blue')
    ax6.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg7['nuv_mag']-sg7['phot_g_mean_mag'], sg7['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax7)
    ax7.axhline(y=horlinelim[i], color='blue')
    ax7.axvline(x=verlinelim[i], color='blue')
    ax7.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    scatter_contour(sg8['nuv_mag']-sg8['phot_g_mean_mag'], sg8['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax8)
    ax8.axhline(y=horlinelim[i], color='blue')
    ax8.axvline(x=verlinelim[i], color='blue')
    ax8.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red')

    ax1.set_xlim((-4, 11))
    ax1.set_ylim((8, -6))
    ax2.set_xlim((-4, 11))
    ax2.set_ylim((8, -6))
    ax3.set_xlim((-4, 11))
    ax3.set_ylim((8, -6))
    ax4.set_xlim((-4, 11))
    ax4.set_ylim((8, -6))
    ax5.set_xlim((-4, 11))
    ax5.set_ylim((8, -6))
    ax6.set_xlim((-4, 11))
    ax6.set_ylim((8, -6))
    ax7.set_xlim((-4, 11))
    ax7.set_ylim((8, -6))
    ax8.set_xlim((-4, 11))
    ax8.set_ylim((8, -6))
    ax1.annotate('gl = 0-45', xy=(-3.5, 7))
    ax2.annotate('gl = 45-90', xy=(-3.5, 7))
    ax3.annotate('gl = 90-135', xy=(-3.5, 7))
    ax4.annotate('gl = 135-180', xy=(-3.5, 7))
    ax5.annotate('gl = 180-225', xy=(-3.5, 7))
    ax6.annotate('gl = 225-270', xy=(-3.5, 7))
    ax7.annotate('gl = 270-315', xy=(-3.5, 7))
    ax8.annotate('gl = 315-360', xy=(-3.5, 7))
    fig.text(0.5, 0.04, 'NUV - G', ha='center')
    fig.text(0.04, 0.5, 'MG', va='center', rotation='vertical')
    plt.suptitle('SEx+G matches, no ext, Hogg par, '+distrange[i]+', iso t=1e9yr, Z=0.0152')
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig('11-29-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_noext.png')
    plt.clf()
    #plt.show()


#######################################################
# HR diagram with extinction by dist and 45 deg slice, colored by ebv
#######################################################
sgcat = fits.open('sex_gaia_dust_interp.fits')[1].data
sga = sgcat[np.where((sgcat['dist'] > 0) & (sgcat['dist'] < 100) & (sgcat['ebv'] > 0))]
sgb = sgcat[np.where((sgcat['dist'] > 100) & (sgcat['dist'] < 300) & (sgcat['ebv'] > 0))]
sgc = sgcat[np.where((sgcat['dist'] > 300) & (sgcat['dist'] < 600) & (sgcat['ebv'] > 0))]
sgd = sgcat[np.where((sgcat['dist'] > 600) & (sgcat['dist'] < 1000) & (sgcat['ebv'] > 0))]
sge = sgcat[np.where((sgcat['dist'] > 1000) & (sgcat['dist'] < 3000) & (sgcat['ebv'] > 0))]
sgf = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['ebv'] > 0))]
cmd = Table.read('cmdfiles/cmd_merged.txt', format='ascii')

sg = [sga, sgb, sgc, sgd, sge]
distrange = ['0 < D < 100 pc', '100 < D < 300 pc', '300 < D < 600 pc', '600 < D < 1000 pc', 'D > 1000 pc']
distfilename = ['0-100', '100-300', '300-600', '600-1000', '1000']
horlinelim = [6.5, 6, 4.5, 3.5, 2.5]
verlinelim = [5, 4, 3.3, 3, 3]


for i in range(len(sg)):
    # By gl
    sg1 = sg[i][np.where((sg[i]['gl_sex'] > 0) & (sg[i]['gl_sex'] < 45))]
    sg2 = sg[i][np.where((sg[i]['gl_sex'] > 45) & (sg[i]['gl_sex'] < 90))]
    sg3 = sg[i][np.where((sg[i]['gl_sex'] > 90) & (sg[i]['gl_sex'] < 135))]
    sg4 = sg[i][np.where((sg[i]['gl_sex'] > 135) & (sg[i]['gl_sex'] < 180))]
    sg5 = sg[i][np.where((sg[i]['gl_sex'] > 180) & (sg[i]['gl_sex'] < 225))]
    sg6 = sg[i][np.where((sg[i]['gl_sex'] > 225) & (sg[i]['gl_sex'] < 270))]
    sg7 = sg[i][np.where((sg[i]['gl_sex'] > 270) & (sg[i]['gl_sex'] < 315))]
    sg8 = sg[i][np.where((sg[i]['gl_sex'] > 315) & (sg[i]['gl_sex'] < 360))]


    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True)
    im1 = ax1.scatter((sg1['nuv_mag']-sg1['ebv']*7.76)-(sg1['phot_g_mean_mag']-sg1['ebv']*3.303), sg1['Mg']-sg1['ebv']*3.303, edgecolor='none', c=np.log10(sg1['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax1.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax1.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax1.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax2.scatter((sg2['nuv_mag']-sg2['ebv']*7.76)-(sg2['phot_g_mean_mag']-sg2['ebv']*3.303), sg2['Mg']-sg2['ebv']*3.303, edgecolor='none', c=np.log10(sg2['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax2.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax2.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax2.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax3.scatter((sg3['nuv_mag']-sg3['ebv']*7.76)-(sg3['phot_g_mean_mag']-sg3['ebv']*3.303), sg3['Mg']-sg3['ebv']*3.303, edgecolor='none', c=np.log10(sg3['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax3.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax3.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax3.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax4.scatter((sg4['nuv_mag']-sg4['ebv']*7.76)-(sg4['phot_g_mean_mag']-sg4['ebv']*3.303), sg4['Mg']-sg4['ebv']*3.303, edgecolor='none', c=np.log10(sg4['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax4.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax4.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax4.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax5.scatter((sg5['nuv_mag']-sg5['ebv']*7.76)-(sg5['phot_g_mean_mag']-sg5['ebv']*3.303), sg5['Mg']-sg5['ebv']*3.303, edgecolor='none', c=np.log10(sg5['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax5.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax5.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax5.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax6.scatter((sg6['nuv_mag']-sg6['ebv']*7.76)-(sg6['phot_g_mean_mag']-sg6['ebv']*3.303), sg6['Mg']-sg6['ebv']*3.303, edgecolor='none', c=np.log10(sg6['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax6.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax6.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax6.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax7.scatter((sg7['nuv_mag']-sg7['ebv']*7.76)-(sg7['phot_g_mean_mag']-sg7['ebv']*3.303), sg7['Mg']-sg7['ebv']*3.303, edgecolor='none', c=np.log10(sg7['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax7.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax7.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax7.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax8.scatter((sg8['nuv_mag']-sg8['ebv']*7.76)-(sg8['phot_g_mean_mag']-sg8['ebv']*3.303), sg8['Mg']-sg8['ebv']*3.303, edgecolor='none', c=np.log10(sg8['ebv']), alpha=0.5, vmin=-1, vmax=0)
    ax8.axhline(y=horlinelim[i], color='black', linewidth=2)
    ax8.axvline(x=verlinelim[i], color='black', linewidth=2)
    ax8.plot(cmd['NUV']-cmd['G'], cmd['G'], c='red', linewidth=2)

    ax1.set_xlim((-4, 11))
    ax1.set_ylim((8, -6))
    ax2.set_xlim((-4, 11))
    ax2.set_ylim((8, -6))
    ax3.set_xlim((-4, 11))
    ax3.set_ylim((8, -6))
    ax4.set_xlim((-4, 11))
    ax4.set_ylim((8, -6))
    ax5.set_xlim((-4, 11))
    ax5.set_ylim((8, -6))
    ax6.set_xlim((-4, 11))
    ax6.set_ylim((8, -6))
    ax7.set_xlim((-4, 11))
    ax7.set_ylim((8, -6))
    ax8.set_xlim((-4, 11))
    ax8.set_ylim((8, -6))
    ax1.annotate('gl = 0-45', xy=(-3.5, 7))
    ax2.annotate('gl = 45-90', xy=(-3.5, 7))
    ax3.annotate('gl = 90-135', xy=(-3.5, 7))
    ax4.annotate('gl = 135-180', xy=(-3.5, 7))
    ax5.annotate('gl = 180-225', xy=(-3.5, 7))
    ax6.annotate('gl = 225-270', xy=(-3.5, 7))
    ax7.annotate('gl = 270-315', xy=(-3.5, 7))
    ax8.annotate('gl = 315-360', xy=(-3.5, 7))
    fig.text(0.5, 0.04, '(NUV - E$_{B-V}$ * 7.76) - (G - E$_{B-V}$ * 3.303)', ha='center')
    fig.text(0.04, 0.5, 'MG - E$_{B-V}$ * 3.303', va='center', rotation='vertical')
    plt.suptitle('SEx+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', iso t=1e9yr, Z=0.0152')
    fig.subplots_adjust(hspace=0, wspace=0)
    fig.subplots_adjust(right=0.9)
    fig.colorbar(im1, cax=fig.add_axes([0.9, 0.15, 0.01, 0.7])).set_label('log10(E$_{B-V}$)')
    plt.savefig('11-29-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'_ebvcolor.png')
    plt.clf()
    #plt.show()



cols = fits.ColDefs([fits.Column(name='tilenum', format='D', array=g['tilenum']), fits.Column(name='nuv_mag', format='D', array=g['nuv_mag']), fits.Column(name='nuv_magerr', format='D', array=g['nuv_magerr']), fits.Column(name='fuv_mag', format='D', array=g['fuv_mag']), fits.Column(name='fuv_magerr', format='D', array=g['fuv_magerr']), fits.Column(name='glon', format='D', array=g['glon']), fits.Column(name='glat', format='D', array=g['glat']), fits.Column(name='ra', format='D', array=g['ra']), fits.Column(name='dec', format='D', array=g['dec'])])
endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('GAISobjall2.fits')




gais = fits.open('gais_total.fits')[1].data
gaia = fits.open('gaia-tycho.fits')[1].data
gaisgal = SkyCoord(gais['ra_gais']*u.deg, gais['dec_gais']*u.deg, frame='icrs')
gaiagal = SkyCoord(gaia['ra']*u.deg, gaia['dec']*u.deg, frame='icrs')
gaisind, gaiaind, angsep, ang3d = search_around_sky(gaisgal, gaiagal, 3*u.arcsec)
comb = hstack([Table(gais[gaisind]), Table(gaia[gaiaind])])

comb.remove_columns(('solution_id', 'source_id', 'random_index', 'ref_epoch', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_good_obs_ac', 'astrometric_n_bad_obs_al', 'astrometric_n_bad_obs_ac', 'astrometric_delta_q', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_primary_flag', 'astrometric_relegation_factor', 'astrometric_weight_al', 'astrometric_weight_ac', 'astrometric_priors_used', 'matched_observations', 'duplicated_source', 'scan_direction_strength_k1', 'scan_direction_strength_k2', 'scan_direction_strength_k3', 'scan_direction_strength_k4', 'scan_direction_mean_k1', 'scan_direction_mean_k2', 'scan_direction_mean_k3', 'scan_direction_mean_k4', 'phot_g_n_obs'))

comb.rename_column('ra', 'ra_tgas')
comb.rename_column('dec', 'dec_tgas')
comb.rename_column('b', 'gb_tgas')
comb.rename_column('l', 'gl_tgas')
comb['angsep'] = angsep

ascii.write(comb, 'gais_tgas.txt', format='ipac')


######################################################################
# CMD with lim mag nuv = 20
######################################################################
mnuv = 20
dist = np.arange(1, 5000, 500)
M = mnuv - 5*np.log10(dist) + 5
g = np.linspace(20, 0, 10)

sg = fits.open('gais_tgas_match_dust.fits')[1].data
sggal = SkyCoord(sg['gl_gais']*u.deg, sg['gb_gais']*u.deg, frame='galactic')

apo = fits.open('APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('BovyRC_TGAS.fits')[1].data
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

scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['MNUV'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))

plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['MNUV'], edgecolor='none',c=c2['dist'], s=40, vmin=0, vmax=2000, label='RC stars')

plt.scatter(c1['nuv_mag']-c1['phot_g_mean_mag'], c1['MNUV'], edgecolor='none', c=c1['dist'], s=40, vmin=0, vmax=2000)

colors = ['']
'''
for i in range(len(dist)):
    plt.axhline(M[i], color='black')
    plt.annotate('D = '+str(dist[i])+' pc', xy=(12, M[i]), size=10)
'''

#plt.scatter(mnuv - g, M,c=dist, vmin=0, vmax=5000, s=80, label='m$_{NUV}$ = 20', marker='s')
plt.legend(scatterpoints=1, loc=3)
plt.xlim((2, 14))
plt.ylim((16, 3))
plt.xlabel('NUV - G')
plt.ylabel('M$_{NUV}$')

plt.colorbar().set_label('dist [pc]')
plt.show()

######################################################################
# CMD with lim mag G = 20
######################################################################
mg = 20
dist = np.arange(1, 5000, 500)
M = mg - 5*np.log10(dist) + 5
g = np.linspace(20, 0, 10)

sg = fits.open('gais_tgas_match_dust.fits')[1].data
sg = sg[~np.isnan(sg['ebv'])]
sggal = SkyCoord(sg['gl_gais']*u.deg, sg['gb_gais']*u.deg, frame='galactic')

apo = fits.open('APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('BovyRC_TGAS.fits')[1].data
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

scatter_contour(sg['nuv_mag']-sg['phot_g_mean_mag'], sg['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['Mg'], edgecolor='none',c=c2['dist'], s=40, vmin=0, vmax=2000, label='RC stars')
plt.scatter(c1['nuv_mag']-c1['phot_g_mean_mag'], c1['Mg'], edgecolor='none', c=c1['dist'], s=40, vmin=0, vmax=2000)

colors = ['']
for i in range(len(dist)):
    plt.axhline(M[i], color='black')
    plt.annotate('D = '+str(dist[i])+' pc', xy=(12, M[i]), size=10)

#plt.scatter(mnuv - g, M,c=dist, vmin=0, vmax=5000, s=80, label='m$_{NUV}$ = 20', marker='s')
plt.legend(scatterpoints=1, loc=3)
plt.xlim((1, 14))
plt.ylim((8, -3))
plt.xlabel('NUV - G')
plt.ylabel('M$_{G}$')
plt.colorbar().set_label('dist [pc]')
plt.show()


######################################################################
# Combine and format RC table
######################################################################
sg = fits.open('gais_tgas_apass_dust.fits')[1].data
sg = sg[~np.isnan(sg['ebv'])]
sggal = SkyCoord(sg['gl_gais']*u.deg, sg['gb_gais']*u.deg, frame='galactic')

apo = fits.open('APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('BovyRC_TGAS.fits')[1].data
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
c1['catalog'] = 'apokask'
c2['catalog'] = 'bovy14'

c1['ALPHAFE'] = c1['ALPHA_M'] * c1['M_H'] / c1['FE_H']

c2['AL_FE'] = c2['AL_H'] / c2['FE_H']
c2['CA_FE'] = c2['CA_H'] / c2['FE_H']
c2['C_FE'] = c2['C_H'] / c2['FE_H']
c2['K_FE'] = c2['K_H'] / c2['FE_H']
c2['MG_FE'] = c2['MG_H'] / c2['FE_H']
c2['MN_FE'] = c2['MN_H'] / c2['FE_H']
c2['NA_FE'] = c2['NA_H'] / c2['FE_H']
c2['NI_FE'] = c2['NI_H'] / c2['FE_H']
c2['N_FE'] = c2['N_H'] / c2['FE_H']
c2['O_FE'] = c2['O_H'] / c2['FE_H']
c2['SI_FE'] = c2['SI_H'] / c2['FE_H']
c2['S_FE'] = c2['S_H'] / c2['FE_H']
c2['TI_FE'] = c2['TI_H'] / c2['FE_H']
c2['V_FE'] = c2['V_H'] / c2['FE_H']



######################################################################
# CMD with lim mag G = 20
######################################################################
#cbarax = 'Fe_H'
#cbarax = 'lnM'
#cbarax = 'lnAge'
cbarax = 'ALPHAFE'

if cbarax == 'Fe_H':
    vmin = -0.5
    vmax = .35
if cbarax == 'lnM':
    vmin = -0.5
    vmax = 1.25
if cbarax == 'lnAge':
    vmin = -1
    vmax = 3
if cbarax == 'ALPHAFE':
    vmin = -0.05
    vmax = 0.3

afeval = []
yaxval = []

for i in range(6, 12):
    cut = np.where(((c2['nuv_mag']-c2['phot_g_mean_mag']) > i) & ((c2['nuv_mag']-c2['phot_g_mean_mag']) < i+1))
    afeval.append(np.mean(c2[cut]['FE_H']))
    yaxval.append(np.mean(c2[cut][cbarax]))
plt.scatter(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['FE_H'], edgecolor='none', c=c2['ALPHAFE'], s=40, vmin=vmin, vmax=vmax, **{"zorder":100})
plt.errorbar(c2['nuv_mag']-c2['phot_g_mean_mag'], c2['FE_H'], xerr=c2['nuv_magerr'], fmt=None, marker=None, mew=0, **{"zorder":0})
plt.scatter([6.5,7.5,8.5,9.5,10.5, 11.5], afeval, s=100, marker='s', c=yaxval, vmin=vmin, vmax=vmax)

if cbarax == 'Fe_H':
    plt.colorbar().set_label('Fe/H')

if cbarax == 'lnM':
    plt.colorbar().set_label('lnM [Msol]')

if cbarax == 'lnAge':
    plt.colorbar().set_label('lnAge [Stellar age]')

if cbarax == 'ALPHAFE':
    plt.colorbar().set_label('Alpha/Fe')

plt.xlabel('NUV - G')
plt.ylabel('Fe/H')
plt.title('GAIS + TGAS, RC stars')
plt.xlim((4,12))
plt.show()

plt.scatter(((c2['nuv_mag']-c2['ebv']*7.24)-(c2['phot_g_mean_mag']-c2['ebv']*3.303))[thin], c2['FE_H'][thin], s=80, c=c2['ALPHAFE'][thin], label='thin', vmin=-0.05, vmax=0.3)

plt.scatter(((c2['nuv_mag']-c2['ebv']*7.24)-(c2['phot_g_mean_mag']-c2['ebv']*3.303))[thick], c2['FE_H'][thick], c=c2['ALPHAFE'][thick], s=80, marker='s', label='thick', vmin=-.05, vmax=0.3)

plt.xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
plt.ylabel('Fe/H')
plt.colorbar().set_label('Alpha/Fe')
plt.legend(scatterpoints=1, loc=2)
plt.show()


# B - V version
plt.scatter(c2['B_AB']-c2['V_AB'], c2['FE_H'], edgecolor='none', c=c2['ALPHAFE'], s=80, vmin=-0.05, vmax=0.3, **{"zorder":100})
plt.errorbar(c2['B_AB']-c2['V_AB'], c2['FE_H'], xerr=c2['B_ABerr']-c2['V_ABerr'], yerr=c2['FE_H_ERR'], fmt=None, marker=None, mew=0, **{"zorder":0})

#plt.scatter((c2['B_AB']-c2['ebv']*3.626)-(c2['V_AB']-c2['ebv']*2.742), c2['FE_H'], edgecolor='none', c=c2['ALPHAFE'], s=80, vmin=-0.05, vmax=0.3, **{"zorder":100})
#plt.errorbar((c2['B_AB']-c2['ebv']*3.626)-(c2['V_AB']-c2['ebv']*2.742), c2['FE_H'], xerr=c2['B_ABerr']-c2['V_ABerr'], yerr=c2['FE_H_ERR'], fmt=None, marker=None, mew=0, **{"zorder":0})

#plt.xlim((-1, 0.4))
plt.xlabel('B - V')
#plt.xlabel('(B - E$_{B-V}$ * 3.626) - (V - E$_{B-V}$ * 2.742)')
plt.ylabel('Fe/H')
#plt.title('GAIS + Bovy RC stars, no extinction')
plt.colorbar().set_label('Alpha/Fe')
plt.show()


######################################################################
# Combine all GAIS tables
######################################################################

 cols = fits.ColDefs([fits.Column(name='nuv_mag', format='D', array=comb['nuv_mag']), fits.Column(name='nuv_magerr', format='D', array=comb['nuv_magerr']), fits.Column(name='fuv_mag', format='D', array=comb['fuv_mag']), fits.Column(name='fuv_magerr', format='D', array=comb['fuv_magerr']), fits.Column(name='gl_gais', format='D', array=comb['gl_gais']), fits.Column(name='gb_gais', format='D', array=comb['gb_gais']), fits.Column(name='ra_gais', format='D', array=comb['ra_gais']), fits.Column(name='dec_gais', format='D', array=comb['dec_gais'])])
endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('GAISobjall2.fits')


galex = fits.open('GAIS_total.fits')[1].data
galexgal = SkyCoord(galex['gl_gais']*u.deg, galex['gb_gais']*u.deg, frame='galactic')
tgas = fits.open('gaia-tycho.fits')[1].data
tgasgal = SkyCoord(tgas['l']*u.deg, tgas['b']*u.deg, frame='galactic')
galexind, tgasind, angsep, ang3d = search_around_sky(galexgal, tgasgal, 3*u.arcsec)
comb = hstack([Table(galex[galexind]),Table(tgas[tgasind])])
comb.rename_column('ra', 'ra_tgas')
comb.rename_column('dec', 'dec_tgas')
comb.rename_column('b', 'gb_tgas')
comb.rename_column('l', 'gl_tgas')
comb['angsep'] = angsep


cols = fits.ColDefs([fits.Column(name='nuv_mag', format='D', array=comb['nuv_mag']), fits.Column(name='nuv_magerr', format='D', array=comb['nuv_magerr']), fits.Column(name='fuv_mag', format='D', array=comb['fuv_mag']), fits.Column(name='fuv_magerr', format='D', array=comb['fuv_magerr']), fits.Column(name='gl_gais', format='D', array=comb['gl_gais']), fits.Column(name='gb_gais', format='D', array=comb['gb_gais']), fits.Column(name='ra_gais', format='D', array=comb['ra_gais']), fits.Column(name='dec_gais', format='D', array=comb['dec_gais']),fits.Column(name='hip', format='J', array=comb['hip']),fits.Column(name='tycho2_id', format='12A', array=comb['tycho2_id']),fits.Column(name='ra_tgas', format='D', array=comb['ra_tgas']),fits.Column(name='ra_error', format='D', array=comb['ra_error']),fits.Column(name='dec_tgas', format='D', array=comb['dec_tgas']),fits.Column(name='dec_error', format='D', array=comb['dec_error']),fits.Column(name='parallax', format='D', array=comb['parallax']),fits.Column(name='parallax_error', format='D', array=comb['parallax_error']),fits.Column(name='pmra', format='D', array=comb['pmra']),fits.Column(name='pmra_error', format='D', array=comb['pmra_error']),fits.Column(name='pmdec', format='D', array=comb['pmdec']),fits.Column(name='pmdec_error', format='D', array=comb['pmdec_error']),fits.Column(name='ra_dec_corr', format='E', array=comb['ra_dec_corr']),fits.Column(name='ra_parallax_corr', format='E', array=comb['ra_parallax_corr']),fits.Column(name='ra_pmra_corr', format='E', array=comb['ra_pmra_corr']),fits.Column(name='ra_pmdec_corr', format='E', array=comb['ra_pmdec_corr']),fits.Column(name='dec_parallax_corr', format='E', array=comb['dec_parallax_corr']),fits.Column(name='dec_pmra_corr', format='E', array=comb['dec_pmra_corr']),fits.Column(name='dec_pmdec_corr', format='E', array=comb['dec_pmdec_corr']),fits.Column(name='parallax_pmra_corr', format='E', array=comb['parallax_pmra_corr']),fits.Column(name='parallax_pmdec_corr', format='E', array=comb['parallax_pmdec_corr']),fits.Column(name='pmra_pmdec_corr', format='E', array=comb['pmra_pmdec_corr']),fits.Column(name='phot_g_mean_flux', format='D', array=comb['phot_g_mean_flux']),fits.Column(name='phot_g_mean_flux_error', format='D', array=comb['phot_g_mean_flux_error']),fits.Column(name='phot_g_mean_mag', format='D', array=comb['phot_g_mean_mag']),fits.Column(name='phot_variable_flag', format='13A', array=comb['phot_variable_flag']),fits.Column(name='gl_tgas', format='D', array=comb['gl_tgas']),fits.Column(name='gb_tgas', format='D', array=comb['gb_tgas']),fits.Column(name='ecl_lon', format='D', array=comb['ecl_lon']),fits.Column(name='ecl_lat', format='D', array=comb['ecl_lat']), fits.Column(name='angsep', format='D', array=comb['angsep'])])
endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('gais_tgas.fits')
vstack


dust = fits.open('gais_tgas_dust.fits')[1].data
tycho = fits.open('tycho2.fits')[1].data
dustgal = SkyCoord(dust['ra_tgas']*u.deg, dust['dec_tgas']*u.deg, frame='icrs')
tychogal = SkyCoord(tycho['RAJ2000']*u.deg, tycho['DEJ2000']*u.deg, frame='icrs')
dustind, tychoind, angsep, ang3d = search_around_sky(dustgal, tychogal, 3*u.arcsec)
d2 = dust[dustind]
t2 = tycho[tychoind]
comb = hstack([Table(d2), Table(t2)])
comb.remove_columns(('TYC1', 'TYC2', 'TYC3', 'HIP', 'DE_ICRS_', 'RA_ICRS_'))
comb['MNUV'] = comb['nuv_mag'] - comb['distmod']
comb['angsep_tgasty'] = angsep
V = comb['VTmag'] - 0.09 * (comb['BTmag']-comb['VTmag'])
B = comb['VTmag'] + 0.085 * (comb['BTmag']-comb['VTmag'])
Verr = comb['e_VTmag'] - 0.09 * (comb['e_BTmag']-comb['e_VTmag'])
Berr = comb['e_VTmag'] + 0.085 * (comb['e_BTmag']-comb['e_VTmag'])
VAB = V - 0.044
BAB = B - 0.163
BABerr = Berr - 0.163
VABerr = Verr - 0.044
comb['B_AB'] = BAB
comb['B_ABerr'] = BABerr
comb['V_AB'] = VAB
comb['V_ABerr'] = VABerr
Gerr = np.sqrt((-2.5/(np.log(10)*comb['phot_g_mean_flux']))**2 * comb['phot_g_mean_flux_error']**2)
comb['Gerr'] = Gerr
#ascii.write(comb, 'gais_tgas_tycho_dust.txt', format='basic')

cols = fits.ColDefs([fits.Column(name='nuv_mag', format='D', array=comb['nuv_mag']), fits.Column(name='nuv_magerr', format='D', array=comb['nuv_magerr']), fits.Column(name='fuv_mag', format='D', array=comb['fuv_mag']), fits.Column(name='fuv_magerr', format='D', array=comb['fuv_magerr']), fits.Column(name='gl_gais', format='D', array=comb['gl_gais']), fits.Column(name='gb_gais', format='D', array=comb['gb_gais']), fits.Column(name='ra_gais', format='D', array=comb['ra_gais']), fits.Column(name='dec_gais', format='D', array=comb['dec_gais']), fits.Column(name='hip', format='K', array=comb['hip']), fits.Column(name='tycho2_id', format='13A', array=comb['tycho2_id']), fits.Column(name='ra_tgas', format='D', array=comb['ra_tgas']), fits.Column(name='ra_error', format='D', array=comb['ra_error']), fits.Column(name='dec_tgas', format='D', array=comb['dec_tgas']), fits.Column(name='dec_error', format='D', array=comb['dec_error']), fits.Column(name='parallax', format='D', array=comb['parallax']), fits.Column(name='parallax_error', format='D', array=comb['parallax_error']), fits.Column(name='pmra', format='D', array=comb['pmra']), fits.Column(name='pmra_error', format='D', array=comb['pmra_error']), fits.Column(name='pmdec', format='D', array=comb['pmdec']), fits.Column(name='pmdec_error', format='D', array=comb['pmdec_error']), fits.Column(name='ra_dec_corr', format='E', array=comb['ra_dec_corr']), fits.Column(name='ra_parallax_corr', format='D', array=comb['ra_parallax_corr']), fits.Column(name='ra_pmra_corr', format='D', array=comb['ra_pmra_corr']), fits.Column(name='ra_pmdec_corr', format='D', array=comb['ra_pmdec_corr']), fits.Column(name='dec_parallax_corr', format='D', array=comb['dec_parallax_corr']), fits.Column(name='dec_pmra_corr', format='D', array=comb['dec_pmra_corr']), fits.Column(name='dec_pmdec_corr', format='D', array=comb['dec_pmdec_corr']), fits.Column(name='parallax_pmra_corr', format='D', array=comb['parallax_pmra_corr']), fits.Column(name='parallax_pmdec_corr', format='D', array=comb['parallax_pmdec_corr']), fits.Column(name='pmra_pmdec_corr', format='D', array=comb['pmra_pmdec_corr']), fits.Column(name='phot_g_mean_mag', format='D', array=comb['phot_g_mean_mag']), fits.Column(name='Gerr', format='D', array=comb['Gerr']), fits.Column(name='phot_g_mean_flux', format='D', array=comb['phot_g_mean_flux']), fits.Column(name='phot_g_mean_flux_error', format='D', array=comb['phot_g_mean_flux_error']), fits.Column(name='phot_variable_flag', format='19A', array=comb['phot_variable_flag']), fits.Column(name='gl_tgas', format='D', array=comb['gl_tgas']), fits.Column(name='gb_tgas', format='D', array=comb['gb_tgas']), fits.Column(name='ecl_lon', format='D', array=comb['ecl_lon']), fits.Column(name='ecl_lat', format='D', array=comb['ecl_lat']), fits.Column(name='angsep', format='D', array=comb['angsep']), fits.Column(name='dist', format='D', array=comb['dist']), fits.Column(name='distmod', format='D', array=comb['distmod']), fits.Column(name='Mg', format='D', array=comb['Mg']), fits.Column(name='ebv', format='D', array=comb['ebv']), fits.Column(name='parallax_hogg', format='D', array=comb['parallax_hogg']), fits.Column(name='B_AB', format='D', array=comb['B_AB']), fits.Column(name='B_ABerr', format='D', array=comb['B_ABerr']), fits.Column(name='V_AB', format='D', array=comb['V_AB']), fits.Column(name='V_ABerr', format='D', array=comb['V_ABerr']), fits.Column(name='MNUV', format='D', array=comb['MNUV']), fits.Column(name='MB', format='D', array=comb['B_AB']-comb['distmod']), fits.Column(name='MV', format='D', array=comb['V_AB']-comb['distmod'])])



##########################################
# Alpha/Fe vs Fe/H vs NUV - G
##########################################
m = (0.095-.21)/(0+0.8)
b = 0.095
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))
color = ((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))

afe_apo = np.where(rc['ALPHA_M'] > 0)
rc['ALPHAFE'][afe_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FE_H'])[afe_apo]
alphafe = rc['ALPHAFE']

# for now use this. err is only for apo
afeerr = np.zeros(len(rc))
afeerr[afe_apo] = (alphafe * np.sqrt((rc['ALPHA_M_err']/rc['ALPHA_M'])**2 + (rc['M_H_err']/rc['M_H'])**2 + (rc['FE_H_err']/rc['FE_H'])**2))[afe_apo]


plt.scatter(rc['FE_H'][thin], rc['ALPHAFE'][thin], s=80, edgecolor='none', c=color[thin], label='Thin disk', vmin=7, vmax=11, marker='D', cmap=cm.jet)

plt.scatter(rc['FE_H'][thick], rc['ALPHAFE'][thick], c=color[thick], s=200, marker='s', label='Thick disk', vmin=7, vmax=11, edgecolor='black', linewidth=3, cmap=cm.jet)

plt.errorbar(rc['FE_H'], rc['ALPHAFE'], xerr=(rc['FE_H_err']), yerr=afeerr, ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})



#plt.xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha$/Fe]')
plt.xlim((-0.7, 0.5))
plt.ylim((-0.05, 0.3))
plt.colorbar().set_label('(NUV - G)$_0$')
leg = plt.legend(scatterpoints=1)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
plt.show()


######################################################################
# Fe/H sv NUV - G with dust correction, errorbars, fit line
######################################################################
from scipy import stats

rc = Table.read('rcall_match.txt', format='ascii')

x1 = rc['nuv_mag'] - rc['phot_g_mean_mag']
x2 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
y = rc['FE_H']

m1, b1, rval1, pval1, stderr1 = stats.linregress(x1, y)
line1 = m1*x1 + b1
err1 = np.sqrt(np.sum((line1-y)**2/len(y)))
m2, b2, rval2, pval2, stderr2 = stats.linregress(x2, y)
line2 = m2*x2 + b2
err2 = np.sqrt(np.sum((line2-y)**2/len(y)))

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]

import scipy.odr as odr
linear = odr.Model(f)
x1odr = odr.Data(x1, y)
x1odr = odr.ODR(x1odr, linear, beta0=[m1, b1])
x1output = x1odr.run()
x2odr = odr.Data(x2, y)
x2odr = odr.ODR(x2odr, linear, beta0=[estm2, estb2])
x2output = x2odr.run()

estm2 = (0.402244+0.485577)/(9.99539-6.59447)


# Separating thin vs thick disk
m = (0.095-.21)/(0+0.8)
b = 0.095
thick, = np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin, = np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
cmap = ax1.scatter((rc['nuv_mag']-rc['phot_g_mean_mag'])[thin], rc['FE_H'][thin], c=rc['ALPHAFE'][thin], s=120, vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, **{"zorder":5})
ax1.errorbar((rc['nuv_mag']-rc['phot_g_mean_mag'])[thin], rc['FE_H'][thin], xerr=(rc['nuv_magerr']-rc['Gerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax1.scatter((rc['nuv_mag']-rc['phot_g_mean_mag'])[thick], rc['FE_H'][thick], c=rc['ALPHAFE'][thick], s=200, vmin=-0.05, vmax=0.3, marker='s', linewidth=3, cmap=plasma,**{"zorder":5})
ax1.errorbar((rc['nuv_mag']-rc['phot_g_mean_mag'])[thick], rc['FE_H'][thick], xerr=(rc['nuv_magerr']-rc['Gerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
#ax1.plot(x1, line1, linewidth=2, c='black', zorder=10)

ax2.scatter(((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thin], rc['FE_H'][thin], c=rc['ALPHAFE'][thin], s=120, label='Thin disk', vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, **{"zorder":5})
ax2.errorbar(((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thin], rc['FE_H'][thin], xerr=(rc['nuv_magerr']-rc['Gerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax2.scatter(((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thick], rc['FE_H'][thick], c=rc['ALPHAFE'][thick], s=200, label='Thick disk', vmin=-0.05, vmax=0.3, marker='s', linewidth=3, cmap=plasma, **{"zorder":5})
ax2.errorbar(((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thick], rc['FE_H'][thick], xerr=(rc['nuv_magerr']-rc['Gerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
#ax2.plot(x2, line2, linewidth=2, c='black', zorder=10)

#ax2.plot(x2, 0.22222906*x2 - 1.92779686, linewidth=2, c='blue', zorder=10)


ax1.set_xlim((5,12))
ax2.set_xlim((5,12))
ax1.set_ylim((-0.7, 0.5))
ax2.set_ylim((-0.7, 0.5))
ax1.set_xlabel('NUV - G')
#ax2.set_xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
ax2.set_xlabel('(NUV - G)$_{0}$')
ax1.set_ylabel('[Fe/H]')
leg = ax2.legend(scatterpoints=1, loc=4)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
ax2.set_xticklabels([' ', '6', '7', '8', '9', '10', '11', '12'])
fig.colorbar(cmap, cax=cbar_ax).set_label(r'[$\alpha$/Fe]')
plt.show()

######################################################################
# Color vs TEFF vs Fe/H
######################################################################
m = (0.095-.21)/(0+0.8)
b = 0.095
thick, = np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin, = np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))
y1 = (rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742)
y2 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
x = rc['TEFF']
m1, b1, rval1, pval1, stderr1 = stats.linregress(x, y1)
line1 = m1*x + b1
err1 = np.sqrt(np.sum((line1-y1)**2/len(y1)))
m2, b2, rval2, pval2, stderr2 = stats.linregress(x, y2)
line2 = m2*x + b2
err2 = np.sqrt(np.sum((line2-y2)**2/len(y2)))

fig, (ax1, ax2) = plt.subplots(2, 1)
cmap = ax1.scatter(rc['TEFF'][thin], ((rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742))[thin], s=80, c=rc['FE_H'][thin], cmap=viridis, vmin=-.5, vmax=.35, marker='D', label='Thin disk')
ax1.scatter(rc['TEFF'][thick], ((rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742))[thick], s=200, c=rc['FE_H'][thick], cmap=viridis, vmin=-.5, vmax=.35, marker='s', label='Thick disk', edgecolor='black', linewidth=3)
ax1.errorbar(rc['TEFF'][thin], ((rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742))[thin], xerr=rc['TEFF_ERR'][thin], yerr=(rc['Berr_apass']-rc['Verr_apass'])[thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax1.plot(x, line1, linewidth=2, c='black', zorder=10)

ax2.scatter(rc['TEFF'][thin], ((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thin], s=80, c=rc['FE_H'][thin], cmap=viridis, vmin=-.5, vmax=.35, marker='D')
ax2.scatter(rc['TEFF'][thick], ((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thick], s=200, c=rc['FE_H'][thick], cmap=viridis, vmin=-.5, vmax=.35, marker='s', linewidth=3)
ax2.errorbar(rc['TEFF'][thin], ((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))[thin], xerr=rc['TEFF_ERR'][thin], yerr=(rc['nuv_magerr']-rc['Gerr'])[thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax2.plot(x, line2, linewidth=2, c='black', zorder=10)

#ax1.set_xlabel('T$_{eff}$')
ax1.set_ylabel('(B - V)$_0$')
ax2.set_xlabel('T$_{eff}$')
ax2.set_ylabel('(NUV - G)$_0$')
ax1.set_xlim((4300, 5200))
ax2.set_xlim((4300, 5200))
ax1.set_ylim((0.5, 1.6))
ax2.set_ylim((5, 12))
fig.subplots_adjust(right=.84)
fig.subplots_adjust(hspace=0)
ax1.set_xticklabels([])
leg = ax1.legend(scatterpoints=1)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
plt.show()


######################################################################
# CMD with MG and MNUV
######################################################################
sg = fits.open('gais_tgas_apass_dust.fits')[1].data
sg = sg[~np.isnan(sg['ebv'])]

fig, (ax1, ax2) = plt.subplots(2, 1)
scatter_contour((sg['nuv_mag']-sg['ebv']*7.24)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['MNUV']-sg['ebv']*7.24), threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=ax2)
scatter_contour((sg['nuv_mag']-sg['ebv']*7.24)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['MG']-sg['ebv']*3.303), threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray), ax=ax1)
ax1.set_xlim((0, 12))
ax1.set_ylim((8, -2))
ax2.set_xlim((0, 12))
ax2.set_ylim((17.4, 0.4))
#ax1.set_xlabel('(NUV - G)$_0$')
ax2.set_xlabel('(NUV - G)$_0$')
ax2.set_ylabel('M$_{NUV_0}$')
ax1.set_ylabel('M$_{G_0}$')
fig.subplots_adjust(hspace=0)
ax1.set_xticklabels([])
#ax2.set_yticklabels(['', '2', '4', '6', '8', '10', '12', '14', '16'])

#ax1.add_patch(matplotlib.patches.Rectangle((7.2, -1.1),2.9,2.4,edgecolor='red',alpha=0.5))
#ax2.add_patch(matplotlib.patches.Rectangle((7.2, 5.9),4.6,1.8,edgecolor='red',alpha=0.5, angle=45))
#ax1.scatter(x[rccut], (sg['MG']-sg['ebv']*3.303)[rccut], edgecolor='none', alpha=0.01)
#ax2.scatter(x[rccut], (sg['MNUV']-sg['ebv']*7.24)[rccut], edgecolor='none', alpha=0.01)
plt.show()


box,= np.where((nuvg > 7.2) & (nuvg < 7.2+2.9) & (sg['MG'] > -1.1 + 2.4) & (sg['MG'] < -1.1))


cols = fits.ColDefs([fits.Column(name='nuv_mag',format='D', array=comb['nuv_mag']),fits.Column(name='gl_gps',format='D', array=comb['gl_gps']),fits.Column(name='gb_gps',format='D', array=comb['gb_gps']),fits.Column(name='ra_gps',format='D', array=comb['ra_gps']),fits.Column(name='dec_gps',format='D', array=comb['dec_gps']),fits.Column(name='ra_tgas',format='D', array=comb['ra_tgas']),fits.Column(name='dec_tgas',format='D', array=comb['dec_tgas']),fits.Column(name='parallax',format='D', array=comb['parallax']),fits.Column(name='parallax_error',format='D', array=comb['parallax_error']),fits.Column(name='phot_g_mean_mag',format='D', array=comb['phot_g_mean_mag']),fits.Column(name='dist',format='D', array=comb['dist']),fits.Column(name='distmod',format='D', array=comb['distmod']),fits.Column(name='MG',format='D', array=comb['MG']),fits.Column(name='ebv',format='D', array=comb['ebv']),fits.Column(name='parallax_hogg',format='D', array=comb['parallax_hogg']),fits.Column(name='Bmag',format='D', array=comb['Bmag']),fits.Column(name='Vmag',format='D', array=comb['Vmag']),fits.Column(name='gmag',format='D', array=comb['gmag']),fits.Column(name='rmag',format='D', array=comb['rmag']),fits.Column(name='imag',format='D', array=comb['imag'])])


endtable = fits.BinTableHDU.from_columns(cols)
endtable.writeto('galexplane_tgas_dust_apass.fits')


avg = []

for i in range(12,21, 1):
    q = np.where((comb['nuv_mag'] > i) & (comb['nuv_mag'] < i+1))
    c2 = comb[q]
    avg.append(np.median(c2['nuv']-c2['nuv_mag']))

#################################################
# Combine all starcat files
#################################################
scans = ['0005', '0014', '0023', '0032', '0041', '0050', '0059', '0068', '0086', '0095', '0104', '0113', '0122', '0140', '0149', '0158', '0167', '0176', '0185', '0194', '0203', '0212', '0221', '0230', '0239', '0248', '0257', '0284', '0293', '0302', '0311', '0320', '0329', '0338', '0347', '0356', '0392', '0428', '0437', '0446', '0455', '0464', '0473', '0482', '0491', '0500', '0671', '0689', '0716', '0743', '0752', '0761', '0770', '0779', '0788', '0797', '0806', '0815', '0824', '0833', '0878', '0887', '0896', '0905', '0914', '0923', '0932', '0941', '0950', '0959', '0968', '0977', '0986', '0995', '1004', '1013', '1022', '1031', '1040', '1049', '1058', '1067', '1076', '1103', '1112', '1121', '1130', '1139', '1148', '1157', '1166', '1193', '1211', '1229', '1238', '1247', '1256', '1265', '1274', '1283', '1292', '1301', '1310', '1319', '1328', '1337', '1346', '1355', '1364', '1373', '1382', '1391', '1400', '1409', '1418', '1427', '1436', '1445', '1454', '1463', '1472', '1481', '1490', '1499', '1508', '1517', '1526', '1535', '1544', '1553', '1562', '1571', '1580', '1589', '1598', '1607', '1616', '1634', '1652', '1661', '1670', '1679', '1724', '1733', '1742', '1751', '1760', '1769', '1778', '1787', '1796', '1805', '1832', '1850', '1904', '1913', '1976', '1985', '2003', '2012', '2030', '2039', '2057', '2066', '2075', '2084', '2093', '2102', '2111', '2120', '2129', '2138', '2147', '2156', '2165', '2174', '2183', '2192', '2201', '2210', '2219', '2228', '2237', '2246', '2255', '2264', '2282', '2291', '2300', '2309', '2318', '2345', '2354', '2363', '2372', '2381', '2390', '2399', '2408', '2417', '2426', '2435', '2444', '2453', '2462', '2471', '2480', '2489', '2498', '2507', '2516', '2525', '2534', '2543', '2552', '2561', '2570', '2588', '2597', '2606', '2615', '2624', '2633', '2642', '2651', '2660', '2669', '2687', '2696', '2705', '2714', '2723', '2732', '2741', '2750', '2759', '2768', '2786', '2795', '2813', '2831', '2840', '2849', '2858', '2867', '2885', '2894', '2903', '2912', '2921', '2930', '2939', '2957', '2975', '2984', '3011', '3020', '3029', '3038', '3047', '3056', '3065', '3083', '3092', '3101', '3155', '3164', '3173', '3182', '3191', '3200', '3209', '3218', '3227', '3236', '3245', '3254', '3263', '3272', '3281', '3290', '3299', '3317', '3326', '3335', '3344', '3353', '3371', '3380', '3389', '3398', '3416', '3425', '3434', '3452', '3488', '3497', '3506', '3515', '3524', '3533', '3542', '3551', '3560', '3578', '3587', '3596']

alldata = Table.read('starcat_0005mapweight_fec_fwhm.txt', format='ascii')
#alldata = alldata[np.where((alldata['FWHM_IMAGE'] < 10) & (alldata['FWHM_IMAGE'] > 3.5))]

alldata['gl'][np.where(alldata['gl'] > 350)] = alldata['gl'][np.where(alldata['gl'] > 350)] - 360
alldata = alldata[np.where(alldata['gl'] < 0.5 + 0.33)]


alldata = Table.read('starcat_1103mapweight_fec_fwhm.txt', format='ascii')

for region in scans[84:152]:
    a = Table.read('starcat_'+region+'mapweight_fec_fwhm.txt', format='ascii')
    number = float(region)/10.
    a = a[np.where(a['gl'] < number + 0.31)]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, 'starcat_allscans_09-05-17.txt', format='basic')



galex = fits.open('GALEXPlane2_smohammed.fit')[1].data
galex['Glon'][np.where(galex['Glon'] > 360)] = galex['Glon'][np.where(galex['Glon'] > 360)] - 360
galexgal = SkyCoord(galex['Glon']*u.deg, galex['Glat']*u.deg, frame='galactic')
catgal = SkyCoord(alldata['gl']*u.deg,alldata['gb']*u.deg, frame='galactic')
catind, galexind, angsep, ang3d = search_around_sky(catgal, galexgal, 3*u.arcsec)
comb = hstack([Table(alldata[catind]), Table(galex[galexind])])
comb['angsep'] = angsep



pscats = ['0-100', '100-200', '200-300', '300-360']

for scan in pscats:
    ps = fits.open('ps1/ps1_planecut_'+scan+'_g10-20.fits')[1].data
    psgal = SkyCoord(ps['glmean']*u.deg, ps['gbmean']*u.deg, frame='galactic')
    catgal = SkyCoord(alldata['gl']*u.deg,alldata['gb']*u.deg, frame='galactic')
    catind, psind, angsep, ang3d = search_around_sky(catgal, psgal, 3*u.arcsec)
    comb = hstack([Table(alldata[catind]), Table(ps[psind])])
    comb['angsep'] = angsep
    ascii.write('starcat_ps1'+str(scan)+'_g10-20.txt', format='basic')



#################################################
# rc temp analysis
#################################################
rc = fits.open('rc_all_10-31.fits')[1].data

nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
q = np.where(x1 < 6.6)
nuvv = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['V_apass']-rc['ebv']*2.742)
nuvb = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['B_apass']-rc['ebv']*3.626)
gv = (rc['phot_g_mean_mag']-rc['ebv']*3.303)-(rc['V_apass']-rc['ebv']*2.742)
gb = (rc['phot_g_mean_mag']-rc['ebv']*3.303)-(rc['B_apass']-rc['ebv']*3.626)
plt.scatter(rc['TEFF'], gb)
plt.scatter(rc['TEFF'][q], gb[q])
plt.show()

bv = (rc['B_apass']-rc['ebv']*3.626)-(rc['v_apass']-rc['ebv']*2.742)


#################################################
# rc fit lines for Fe/H vs NUV - G vs Alpha/Fe
#################################################
rc = fits.open('rc_all_10-31.fits')[1].data
xa1 = rc['nuv_mag']-rc['phot_g_mean_mag']
xb1 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
y1 = rc['FE_H']

q = np.where(xb1 > 6.6)
w = np.where(xb1 < 6.6)
xa2 = xa1[q]
xb2 = xb1[q]
y2 = y1[q]

za1, va1 = np.polyfit(xa1, y1, 1, cov=True)
zb1, vb1 = np.polyfit(xb1, y1, 1, cov=True)
za2, va2 = np.polyfit(xa2, y2, 1, cov=True)
zb2, vb2 = np.polyfit(xb2, y2, 1, cov=True)

pa1 = np.poly1d(za1)
pb1 = np.poly1d(zb1)
pa2 = np.poly1d(za2)
pb2 = np.poly1d(zb2)
xp = np.linspace(6, 11, 50)

a1err = np.sum(np.sqrt(np.diag(va1)))
a2err = np.sum(np.sqrt(np.diag(va2)))
b1err = np.sum(np.sqrt(np.diag(vb1)))
b2err = np.sum(np.sqrt(np.diag(vb2)))

'''
# Now get std errors
ma1, ba1, rvala1, pvala1, stderra1 = stats.linregress(xa1, y1)
ma2, ba2, rvala2, pvala2, stderra2 = stats.linregress(xa2, y2)
mb1, bb1, rvalb1, pvalb1, stderrb1 = stats.linregress(xb1, y1)
mb2, bb2, rvalb2, pvalb2, stderrb2 = stats.linregress(xb2, y2)
'''

m = (0.095-.21)/(0+0.8)
b = 0.095
thick, = np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin, = np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
cmap = ax1.scatter(xa1[thin], y1[thin], c=rc['ALPHAFE'][thin], s=120, vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, **{"zorder":5})
ax1.errorbar(xa1[thin], y1[thin], xerr=(rc['nuv_magerr']-rc['Gerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax1.scatter(xa1[thick], y1[thick], c=rc['ALPHAFE'][thick], s=200, vmin=-0.05, vmax=0.3, marker='s', edgecolor='blue', linewidth=2, cmap=plasma,**{"zorder":5})
ax1.errorbar(xa1[thick], y1[thick], xerr=(rc['nuv_magerr']-rc['Gerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax1.plot(xp, pa1(xp), linewidth=2, c='red', zorder=10) # no dust, all points
ax1.plot(xp, pa2(xp), linewidth=2, c='black', zorder=10) # no dust, only black points
ax1.scatter(xa1[w], y1[w], c=rc['ALPHAFE'][w], s=120, vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, edgecolor='red', linewidth=2, **{"zorder":5})

ax2.scatter(xb1[thin], y1[thin], c=rc['ALPHAFE'][thin], s=120, label='Thin disk', vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, **{"zorder":5})
ax2.errorbar(xb1[thin], y1[thin], xerr=(rc['nuv_magerr']-rc['Gerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax2.scatter(xb1[thick], y1[thick], c=rc['ALPHAFE'][thick], s=200, label='Thick disk', vmin=-0.05, vmax=0.3, marker='s', edgecolor='blue', linewidth=2, cmap=plasma, **{"zorder":5})
ax2.errorbar(xb1[thick], y1[thick], xerr=(rc['nuv_magerr']-rc['Gerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax2.plot(xp, pb1(xp), linewidth=2, c='red', zorder=10)
ax2.plot(xp, pb2(xp), linewidth=2, c='black', zorder=10)
ax2.scatter(xb1[w], y1[w], c=rc['ALPHAFE'][w], s=120, vmin=-0.05, vmax=0.3, marker='D', cmap=plasma, edgecolor='red', linewidth=2, **{"zorder":5})

ax1.set_xlim((5,11.9))
ax2.set_xlim((5,11.9))
ax1.set_ylim((-0.7, 0.5))
ax2.set_ylim((-0.7, 0.5))
ax1.set_xlabel('NUV - G')
#ax2.set_xlabel('(NUV - E$_{B-V}$ * 7.24) - (G - E$_{B-V}$ * 3.303)')
ax2.set_xlabel('(NUV - G)$_{0}$')
ax1.set_ylabel('[Fe/H]')
leg = ax2.legend(scatterpoints=1, loc=4)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
ax1.annotate('[Fe/H] = 0.172 * (NUV-G) - 1.552', xy=(5.1, 0.43), color='red', size=15)
ax1.annotate('[Fe/H] = 0.180 * (NUV-G) - 1.630 ', xy=(5.1, 0.47), color='black', size=15)
ax2.annotate('[Fe/H] = 0.203 * (NUV-G)$_0$ - 1.751', xy=(5.1, 0.43), color='red', size=15)
ax2.annotate('[Fe/H] = 0.230 * (NUV-G)$_0$ - 1.987', xy=(5.1, 0.47), color='black', size=15)

ax1.annotate('$\sigma$ = 0.112', xy=(5.1, 0.35), color='red', size=15) # a1
ax1.annotate('$\sigma$ = 0.117', xy=(5.1, 0.39), color='black', size=15) # a2
ax2.annotate('$\sigma$ = 0.094', xy=(5.1, 0.35), color='red', size=15) # b1
ax2.annotate('$\sigma$ = 0.096', xy=(5.1, 0.39), color='black', size=15) # b2

fig.colorbar(cmap, cax=cbar_ax).set_label(r'[$\alpha$/Fe]')
plt.show()



#################################################
# photutils parameters
#################################################
import photutils
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

img = fits.open('im1_0014.fits')[0].data
mean, median, std = sigma_clipped_stats(img, sigma=3.0, iters=5)  
daofind = DAOStarFinder(fwhm=3.0, threshold=3.*std)
sources = daofind(img - median)



from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = photutils.CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch())
plt.imshow(img, cmap=cm.gray, origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)

from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

sigma_clip = SigmaClip(sigma=3., iters=10)
bkg_est = MedianBackground()

bkg = Background2D(img, (30, 30), filter_size=(7,7), sigma_clip=sigma_clip, bkg_estimator=bkg_est)



cmd = Table.read('cmdfiles/cmd_merged_zt.txt', format='ascii')
zrange = np.unique(cmd['Z'])
agerange = np.unique(cmd['logage'])

for i in zrange:
	for j in agerange:
		scatter_contour((sg['nuv_mag']-sg['ebv']*7.24)-(sg['phot_g_mean_mag']-sg['ebv']*3.303), (sg['MG']-sg['ebv']*3.303), threshold=1000, log_counts=True, histogram2d_args=dict(bins=40),plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray))
		cmap = plt.scatter((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303), (rc['MG']-rc['ebv']*3.303), c=rc['FE_H'], vmin=-0.5, vmax=0.35)
		plt.colorbar(cmap).set_label('[Fe/H]')
		plt.xlim(0, 12)
		plt.ylim(8, -2)
		cmd1 = cmd[np.where((cmd['Z'] == i) & (cmd['logage'] == j))]
		plt.scatter(cmd1['NUV']-cmd1['G'], cmd1['G'], c='red', s=1)
		plt.xlabel('(NUV - G)$_{0}$')
		plt.ylabel('M$_{G_0}$')
		plt.title('Z = '+str(i)+', logage = '+str(j))
		plt.savefig('images/12-08-gaiacmd/Z'+str(i)+'-logage'+str(j)+'.png')
		plt.clf()

################################################################
# Calculate errors, get table for paper
################################################################

rc = fits.open('rc_all_10-31.fits')[1].data
ra = rc['ra_tgas']
dec = rc['dec_tgas']
nuv = rc['nuv_mag'] - rc['ebv']*7.24
nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
nuverr = rc['nuv_magerr']
nuvgerr = np.sqrt(rc['nuv_magerr']**2-rc['Gerr']**2)
B = rc['B_apass'] - rc['ebv']*3.626
bv = (rc['B_apass']-rc['ebv']*3.626)-(rc['v_apass']-rc['ebv']*2.742)
Berr = rc['Berr_apass']
bverr = np.sqrt(rc['Berr_apass']**2 - rc['Verr_apass']**2)
bverr[np.isnan(bverr)] = Berr[np.isnan(bverr)]
ebv = rc['ebv']
dm = rc['distmod']
feh = rc['FE_H']
feherr = rc['FE_H_err']
teff = rc['TEFF']
tefferr = rc['TEFF_err']

afe_apo = np.where(rc['ALPHA_M'] > 0)
rc['ALPHAFE'][afe_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FE_H'])[afe_apo]
alphafe = rc['ALPHAFE']

# for now use this. err is only for apo
afeerr = np.zeros(len(rc))
afeerr[afe_apo] = (alphafe * np.sqrt((rc['ALPHA_M_err']/rc['ALPHA_M'])**2 + (rc['M_H_err']/rc['M_H'])**2 + (rc['FE_H_err']/rc['FE_H'])**2))[afe_apo]



ra = ['{:.4f}'.format(x) for x in ra]
dec = ['{:.4f}'.format(x) for x in dec]
nuv = ['{:.2f}'.format(x) for x in nuv]
nuverr = ['{:.2f}'.format(x) for x in nuverr]
nuvg = ['{:.2f}'.format(x) for x in nuvg]
nuvgerr = ['{:.2f}'.format(x) for x in nuvgerr]
B = ['{:.2f}'.format(x) for x in B]
bv = ['{:.2f}'.format(x) for x in bv]
Berr = ['{:.2f}'.format(x) for x in Berr]
bverr = ['{:.2f}'.format(x) for x in bverr]
ebv = ['{:.2f}'.format(x) for x in ebv]
dm = ['{:.2f}'.format(x) for x in dm]
feh = ['{:.2f}'.format(x) for x in feh]
feherr = ['{:.2f}'.format(x) for x in feherr]
teff = ['{:.2f}'.format(x) for x in teff]
alphafe = ['{:.2f}'.format(x) for x in alphafe]
afeerr = ['{:.2f}'.format(x) for x in afeerr]



table = Table([ra, dec, nuv, nuverr, nuvg, nuvgerr, B, Berr, bv, bverr, ebv, dm, feh, feherr, teff, alphafe, afeerr])

t1 = table[:10]
ascii.write(t1, format='latex')


################################################################
# Make plots for scan outputs
################################################################
scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '189-199', '198-208', '207-217', '216-226', '270-280', '279-289', '288-298', '297-307', '351-1', '72-82', '144-154', '225-235', '306-316', '81-91', '153-163', '234-244', '315-325', '90-100', '162-172', '243-253', '324-334', '9-19', '171-181', '252-262', '333-343', '99-109', '180-190', '261-271', '342-352']

for scan in scans:
	t = Table.read('starcat_'+scan+'_01_10_2018.txt', format='ascii')
        t = t[np.where((t['FWHM_IMAGE'] > 0) & (t['FWHM_IMAGE'] < 10))]
	plt.scatter(t['gl'], t['gb'], c=t['FWHM_IMAGE'], vmin=0, vmax=10, s=3)
	plt.colorbar().set_label('FWHM')
	plt.xlabel('gl')
	plt.ylabel('gb')
        plt.title(scan+', 0 > FWHM > 10')
	plt.savefig('01-10-18_glgbvsfwhm_'+scan+'fwhmcut_0-10.png')
	plt.clf()

for scan in scans:
	t = Table.read('starcat_'+scan+'_01_10_2018.txt', format='ascii')
	plt.scatter(t['nuv'], t['FWHM_IMAGE'], s=1, c='black')
	plt.ylim(0, 10)
	plt.savefig('01-10-18_fwhmvsnuv_'+scan+'.png')
	plt.clf()

from scan in scans:
	hdu = fits.open('count_map'+scan+'_in.fits')[0]
	wcs = WCS(hdu.header)
	fig = plt.figure()
	fig.add_subplot(111, projection=wcs)
	plt.imshow(hdu.data, origin='lower', cmap=cm.gray, vmin=0, vmax=0.1)
	plt.title(scan)
	plt.ylim(-11, 11)
	plt.savefig('01-10-18_glgb_'+scan+'.png')
	plt.clf()



scans = ['63-73', '180-190', '225-235']

for scan in scans:
	# Get cmed and imed
	ct = fits.open('im1_'+scan+'_count.fits')[0].data
	ct1 = ct[np.where(ct > 0)]
	img = fits.open('im1_'+scan+'.fits')[0].data
	im1 = img[np.where(img > 10**-5)]
	cmed = np.median(ct1)
	imed = np.median(im1)
	n = imed * (1+(np.random.poisson(cmed*10, img.shape)/10. - cmed)/cmed)

	# make mask and apply to img
	cutoff = np.where(img < 10**-5)
	zeromask = np.zeros(img.shape)
	zeromask[cutoff] = True
	zeromask = zeromask + n

	im2 = img + zeromask

	

# for 63-73
cmed = 0.29931211
imed = 0.039932579


