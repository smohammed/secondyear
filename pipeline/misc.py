####################################################################
# NUV - J vs J - K
####################################################################
star = Table.read('newfield_2mass_t2_jlim_13.5_3arcsec.txt',format='ascii')
newt = Table.read('galex0data_2mass.txt', format='ascii')
pickles = Table.read('picklemags_laphare.txt', format='ascii')

gbrange = 5.

scut = np.where((star['gb_sex'] > -5) & (star['gb_sex'] < 0))
scut2 = np.where((star['gb_sex'] > 0) & (star['gb_sex'] < 5))
ncut = np.where((newt['gb_galex'] > -5) & (newt['gb_galex'] < 0))
ncut2 = np.where((newt['gb_galex'] > 0) & (newt['gb_galex'] < 5))

f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})

ax1.set_title('3", J < 13.5,    VPHAS+T2, -5 < gb < 0')
ax1.set_xlim((-0.5,2.5))
ax1.set_ylim((0,14))
ax1.set_xlabel('J - K')
ax1.set_ylabel('NUV - J')
#ax1.plot([-1.0,1.5],[-1.5,14])
ax1.axhline(y=4,color='black')

a2 = ax1.scatter(newt['j'][ncut]-newt['k'][ncut],newt['nuv_mag'][ncut]-newt['j'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)

a1 = ax1.scatter(star['j'][scut]-star['k'][scut],star['nuv'][scut]-star['j'][scut],edgecolor='none',alpha=0.3)

ax1.legend([a1,a2],['Sextractor','GALEX'],scatterpoints=1,loc=2)


ax2.set_title('3", J < 13.5, 5 < gb < 10')
ax2.set_xlim((-0.5,2.5))
ax2.set_ylim((0,14))
ax2.set_xlabel('J - K')
ax2.set_ylabel('NUV - J')
#ax2.plot([-1.0,1.5],[-1.5,14])
ax2.axhline(y=4,color='black')
#b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
b2 = ax2.scatter(newt['j'][ncut2]-newt['k'][ncut2],newt['nuv_mag'][ncut2]-newt['j'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)

b1 = ax2.scatter(star['j'][scut2]-star['k'][scut2],star['nuv'][scut2]-star['j'][scut2],edgecolor='none',alpha=0.3)
ax2.legend([b1,b2],['Sextractor','GALEX'],scatterpoints=1,loc=2)
plt.show()

####################################################################
# J-H vs H-K
####################################################################
star = Table.read('newfield_2mass_t2_jlim_13.5_3arcsec.txt',format='ascii')
newt = Table.read('galex0data_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags.txt', format='ascii')

gbrange = 5.
scut = np.where((np.abs(star['gb_sex']) > gbrange))
scut2 = np.where((np.abs(star['gb_sex']) < gbrange))
ncut = np.where(np.abs(newt['gb_galex']) > gbrange)
ncut2 = np.where(np.abs(newt['gb_galex']) < gbrange)

f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})

ax1.set_title('3", J < 13.5, Deeper Field, gb > '+str(gbrange))
ax1.set_xlim((-0.4,0.6))
ax1.set_ylim((-0.3,1.0))
ax1.set_xlabel('H - K')
ax1.set_ylabel('J - H')
ax1.plot([-0.4,0.6],[-0.3,1.0])
a2 = ax1.scatter(newt['h'][ncut]-newt['k'][ncut],newt['j'][ncut]-newt['h'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)
a3 = ax1.scatter(pickles['h']-pickles['k'],pickles['j']-pickles['h'],c='black',s=5)
for j in range(len(pickles)):
    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['h'][j]-pickles['k'][j],pickles['j'][j]-pickles['h'][j]),size=12)

a1 = ax1.scatter(star['h'][scut]-star['k'][scut],star['j'][scut]-star['h'][scut],edgecolor='none',alpha=0.1)

ax1.legend([a1,a2,a3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)
ax2.set_title('3", J < 13.5, Deeper Field, gb < '+str(gbrange))
ax2.set_xlim((-0.4,0.6))
ax2.set_ylim((-0.3,1.0))
ax2.set_xlabel('H - K')
ax2.set_ylabel('J - H')
ax2.plot([-0.4,0.6],[-0.3,1.0])
b2 = ax2.scatter(newt['h'][ncut2]-newt['k'][ncut2],newt['j'][ncut2]-newt['h'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)
b3 = ax2.scatter(pickles['h']-pickles['k'],pickles['j']-pickles['h'],c='black',s=5)
for j in range(len(pickles)):
    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['h'][j]-pickles['k'][j],pickles['j'][j]-pickles['h'][j]),size=12)

b1 = ax2.scatter(star['h'][scut2]-star['k'][scut2],star['j'][scut2]-star['h'][scut2],edgecolor='none',alpha=0.1)
ax2.legend([b1,b2,b3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)
plt.show()


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
cm.set_label('Limiting radius [arcsec]')
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

####################################################################
# Compare NUV
####################################################################
# Only detections
sex = Table.read('newfield_gal_ipac_tests.txt',format='ascii')
sexgal = SkyCoord(sex['ra']*u.degree, sex['dec']*u.degree, frame='icrs')
galex = Table.read('galex0data.txt',format='ascii')
galexgal = SkyCoord(galex['ra']*u.deg,galex['dec']*u.deg,frame='icrs')
sexind, galexind, angsep, dist3d = search_around_sky(sexgal,galexgal,5*u.arcsec)
s2 = sex[sexind]
g2 = galex[galexind]

# With 2MASS data
star = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt13.5_tests.txt',format='ascii')
galex2m = Table.read('galex0data_2mass.txt',format='ascii')
stargal = SkyCoord(star['ra_2mass']*u.deg,star['dec_2mass']*u.deg,frame='icrs')
galex2mgal = SkyCoord(galex2m['ra_2mass']*u.deg,galex2m['dec_2mass']*u.deg,frame='icrs')
starind, galex2mind, angsep,dist3d = search_around_sky(stargal,galex2mgal,1*u.arcsec)
star2 = star[starind]
g2mass = galex2m[galex2mind]

plt.scatter(g2['nuv_mag'],s2['nuv']-g2['nuv_mag'],edgecolor='none',alpha=0.5,label='Detections')
plt.scatter(g2mass['nuv_mag'],star2['nuv']-g2mass['nuv_mag'],edgecolor='none',facecolor='red',alpha=0.5,label='2MASS matches')
plt.xlabel('GALEX NUV')
plt.ylabel('SExtractor - GALEX NUV')
plt.xlim((10,24))
plt.ylim((-6,3))
plt.axhline(y=0,c='black')
plt.legend(loc=3,scatterpoints=1)
plt.show()

######################################################################
# GALEX, Pickles, SExtractor with extinction
######################################################################
nuv = ['nuv_mag', 'nuv', 'nuv', 2.9720]
b = ['BJmag', 'b', 'BJmag', 1.3429]
v = ['VJmag', 'v', 'VJmag', 1.0]
j = ['j', 'j', 'j', 0.2876]
h = ['h', 'h', 'h', 0.1783]
k = ['k', 'k', 'k', 0.1170]
u = ['u','u','u',1.5916]
g = ['g','g','g',1.1838]
r = ['r','r','r',0.8664]
ib  = ['i','i','i',0.6418]

plt.scatter(sexv['g_AB']-sexv['i_AB'],sexv['nuv']-sexv['g_AB'],edgecolor='none',facecolor='blue',alpha=0.3,label='VPHAS')
plt.scatter(starv['g_AB']-starv['i_AB'],starv['nuv']-starv['g_AB'],edgecolor='none',facecolor='red',label='VPHAS + 2MASS',alpha=0.3)
plt.xlim((-2,3))
plt.ylim((-2,8))
plt.ylabel('NUV - g (ABmag)')
plt.xlabel('g - i (ABmag)')
plt.legend(loc=3,scatterpoints=1)
plt.gca().invert_yaxis()
plt.scatter(pickles['g']-pickles['i'], pickles['nuv']-pickles['g'], c='black', s=10)

plt.arrow(-1.5, 4,g[3]-ib[3], nuv[3]-g[3], head_length=0.1,head_width=0.07,color='black')

for j in range(0,len(pickles),10):
        plt.annotate(pickles['name'][j], xy=((pickles['g'][j]-pickles['i'][j])+0.01, (pickles['nuv'][j]-pickles['g'][j]) - 0.07), size=15)
plt.show()

delang = 2*np.arcsin(np.sqrt(np.sin((star['dec_sex']-star['dec_2mass'])/2)**2+np.cos(star['dec_sex'])*np.cos(star['dec_2mass'])*np.sin((star['ra_sex']-star['ra_2mass'])/2)**2))

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


labels = ['NUV-i','u-i','g-i','r-i','J-i','H-i','K-i']
plt.xticks([1,2,3,4,5,6,7],labels)
for i in range(len(wd2m)):
    plt.plot([1,2,3,4,5,6,7],cols[[i,i+len(wd2m),i+len(wd2m)*2,i+len(wd2m)*3,i+len(wd2m)*4,i+len(wd2m)*5,i+len(wd2m)*6]],alpha=0.01,color='black')
for i in range(len(wd2m[cut])):
    plt.plot([1,2,3,4,5,6,7],cols2[[i,i+len(wd2m[cut]),i+len(wd2m[cut])*2,i+len(wd2m[cut])*3,i+len(wd2m[cut])*4,i+len(wd2m[cut])*5,i+len(wd2m[cut])*6]],alpha=0.05,color='red')

plt.plot([1,2,3,4,5,6,7],cols[[1,1+len(wd2m),1+len(wd2m)*2,1+len(wd2m)*3,1+len(wd2m)*4,1+len(wd2m)*5,1+len(wd2m)*6]],alpha=0.1,color='black',label='VPHAS')
plt.plot([1,2,3,4,5,6,7],cols2[[1,1+len(wd2m[cut]),1+len(wd2m[cut])*2,1+len(wd2m[cut])*3,1+len(wd2m[cut])*4,1+len(wd2m[cut])*5,1+len(wd2m[cut])*6]],alpha=0.1,color='red',label='VPHAS + 2MASS')
plt.legend(loc=2)
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda$ - i (ABmag)')
plt.ylim((5,-12))
plt.title('VPHAS + 2MASS, Verbeek WDs')
plt.show()

####################################################################
# GALEX and VPHAS sky plot
####################################################################
galex = fits.open('GALEXAIS.fits')[1].data
galexcut = np.where((galex['glon'] > 0) & (galex['glon'] < 40))
g2 = galex[galexcut]
galexgal = SkyCoord(g2['glon']*u.deg,g2['glat']*u.deg,frame='galactic')
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

v1cut = np.where((v1['gl'] > 0) & (v1['gl'] < 40))
v1 = v1[v1cut]

wd = Table.read('wds_gais_vphas_newgrcut.txt',format='ascii')
wdcut = np.where((wd['gl_galex'] > 0) & (wd['gl_galex'] < 40))
wd2 = wd[wdcut]

# SED plot
v1nuvi = np.zeros(len(v1))
v1ui = v1['u_AB'] - v1['i_AB']
v1gi = v1['g_AB'] - v1['i_AB']
v1ri = v1['r_AB'] - v1['i_AB']

nuvi = wd2['nuv_mag'] - wd2['i_AB']
ui = wd2['u_AB'] - wd2['i_AB']
gi = wd2['g_AB'] - wd2['i_AB']
ri = wd2['r_AB'] - wd2['i_AB']

colsv1 = np.ones(len(v1)*4)
colsv1[0+len(v1)*0:len(v1)*1] = v1nuvi
colsv1[0+len(v1)*1:len(v1)*2] = v1ui
colsv1[0+len(v1)*2:len(v1)*3] = v1gi
colsv1[0+len(v1)*3:len(v1)*4] = v1ri

cols = np.ones(len(wd2)*4)
cols[0+len(wd2)*0:len(wd2)*1] = nuvi
cols[0+len(wd2)*1:len(wd2)*2] = ui
cols[0+len(wd2)*2:len(wd2)*3] = gi
cols[0+len(wd2)*3:len(wd2)*4] = ri

labels = ['NUV-i','u-i','g-i','r-i']
plt.xticks([1,2,3,4],labels)
# All VPHAS
for i in range(len(v1)):
    plt.plot([1,2,3,4],colsv1[[i,i+len(v1),i+len(v1)*2,i+len(v1)*3]],alpha=0.05,color='black')

# VPHAS + GAIS
for i in range(len(wd2)):
    plt.plot([1,2,3,4],cols[[i,i+len(wd2),i+len(wd2)*2,i+len(wd2)*3]],alpha=0.1,color='red')

# Plot this to make legend
plt.plot([1,2,3,4],colsv1[[1,1+len(v1),1+len(v1)*2,1+len(v1)*3]],alpha=0.4,color='black',label='VPHAS only')
plt.plot([1,2,3,4],cols[[1,1+len(wd2),1+len(wd2)*2,1+len(wd2)*3]],alpha=0.4,color='red',label='VPHAS + GAIS')
plt.legend(loc=4)
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda$ - i (ABmag)')
plt.ylim((5,-1.5))
plt.title('VPHAS + GAIS WDs, new gr cut, 0 < gl < 40')
plt.show()


######################################################################
# lambda - i histograms
######################################################################
v = Table.read('wds_vphasonly.txt',format='ascii')
wd = Table.read('wds_gais_vphas_newgrcut.txt',format='ascii')
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
axes[0,0].hist((v['r_AB']-v['i_AB'])[v1cut],bins=bins,label='VPHAS only WDs')
axes[0,1].hist((v['r_AB']-v['i_AB'])[v2cut],bins=bins)
axes[1,0].hist((v['r_AB']-v['i_AB'])[v3cut],bins=bins)
axes[1,1].hist((v['r_AB']-v['i_AB'])[v4cut],bins=bins)
axes[0,0].hist((wd['r_AB']-wd['i_AB'])[wd1cut],bins=bins, label='VPHAS+GAIS WDs')
axes[0,1].hist((wd['r_AB']-wd['i_AB'])[wd2cut],bins=bins)
axes[1,0].hist((wd['r_AB']-wd['i_AB'])[wd3cut],bins=bins)
axes[1,1].hist((wd['r_AB']-wd['i_AB'])[wd4cut],bins=bins)
axes[0,0].annotate('0 < gl < 40',xy=(1,400),size=15)
axes[0,1].annotate('200 < gl < 250',xy=(1,400),size=15)
axes[1,0].annotate('250 < gl < 300',xy=(1,400),size=15)
axes[1,1].annotate('300 < gl < 360',xy=(1,400),size=15)
fig.subplots_adjust(wspace=0,hspace=0)
fig.suptitle(r'WDs, $\Delta \theta_{max}$ = 1.5", gr cut')
axes[1,0].set_xlabel('r - i (ABmag)')
axes[1,1].set_xlabel('r - i (ABmag)')
axes[0,0].legend(scatterpoints=1)
plt.show()

# NUV only
bins= np.linspace(-3,5,20)
fig,axes = plt.subplots(2,2,sharey=True,sharex=True)
axes[0,0].hist((wd['nuv_mag']-wd['i_AB'])[wd1cut],bins=bins, label='VPHAS+GAIS WDs')
axes[0,1].hist((wd['nuv_mag']-wd['i_AB'])[wd2cut],bins=bins)
axes[1,0].hist((wd['nuv_mag']-wd['i_AB'])[wd3cut],bins=bins)
axes[1,1].hist((wd['nuv_mag']-wd['i_AB'])[wd4cut],bins=bins)
axes[0,0].annotate('0 < gl < 40',xy=(2,80),size=15)
axes[0,1].annotate('200 < gl < 250',xy=(2,80),size=15)
axes[1,0].annotate('250 < gl < 300',xy=(2,80),size=15)
axes[1,1].annotate('300 < gl < 360',xy=(2,80),size=15)

fig.subplots_adjust(wspace=0,hspace=0)

fig.suptitle(r'WDs, $\Delta \theta_{max}$ = 1.5", gr cut')

axes[1,0].set_xlabel('NUV - i (ABmag)')
axes[1,1].set_xlabel('NUV - i (ABmag)')
axes[0,0].legend(scatterpoints=1)
plt.show()





