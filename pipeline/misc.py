####################################################################
# Made pickles SED dataset
####################################################################
name,nuv,b,v,j,h,k = [],[],[],[],[],[],[]
picklefiles = np.loadtxt('PICKLES/pickles.txt', dtype='string')
for i in picklefiles:
    a = Table.read('PICKLES/'+i,format='ascii')
    name.append(i)
    nuv.append(a['col2'][224])
    b.append(a['col2'][641])
    v.append(a['col2'][781])
    j.append(a['col2'][2241])
    h.append(a['col2'][3095])
    k.append(a['col2'][4089])

####################################################################
# NUV - B vs B - V plot
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
ax1.set_xlim((-0.5,2))
ax1.set_ylim((2,10))
ax1.set_xlabel('B - V')
ax1.set_ylabel('NUV - B')
ax1.plot([-0.5,2],[2,10])
a2 = ax1.scatter(newt['BJmag'][ncut]-newt['VJmag'][ncut],newt['nuv_mag'][ncut]-newt['BJmag'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)
a3 = ax1.scatter(pickles['b']-pickles['v'],pickles['nuv']-pickles['b'],c='black',s=5)
for j in range(len(pickles)):
    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['b'][j]-pickles['v'][j],pickles['nuv'][j]-pickles['b'][j]),size=12)

a1 = ax1.scatter(star['BJmag'][scut]-star['VJmag'][scut],star['nuv'][scut]-star['BJmag'][scut],edgecolor='none',alpha=0.5)
ax1.legend([a1,a2,a3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=4)

ax2.set_title('3", J < 13.5, Deeper Field, gb < '+str(gbrange))
ax2.set_xlim((-0.5,2))
ax2.set_ylim((2,10))
ax2.set_xlabel('B - V')
ax2.set_ylabel('NUV - B')
ax2.plot([-0.5,2],[2,10])
b2 = ax2.scatter(newt['BJmag'][ncut2]-newt['VJmag'][ncut2],newt['nuv_mag'][ncut2]-newt['BJmag'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)
b3 = ax2.scatter(pickles['b']-pickles['v'],pickles['nuv']-pickles['b'],c='black',s=5)
for j in range(len(pickles)):
    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['b'][j]-pickles['v'][j],pickles['nuv'][j]-pickles['b'][j]),size=12)

b1 = ax2.scatter(star['BJmag'][scut2]-star['VJmag'][scut2],star['nuv'][scut2]-star['BJmag'][scut2],edgecolor='none',alpha=0.5)
ax2.legend([b1,b2,b3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=4)
plt.show()

    

####################################################################
# NUV - J vs J - K
####################################################################
star = Table.read('newfield_2mass_t2_jlim_13.5_3arcsec.txt',format='ascii')
newt = Table.read('galex0data_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags_laphare.txt', format='ascii')

gbrange = 5.

scut = np.where((star['gb_sex'] > -5) & (star['gb_sex'] < 0))
scut2 = np.where((star['gb_sex'] > 0) & (star['gb_sex'] < 5))
ncut = np.where((newt['gb_galex'] > -5) & (newt['gb_galex'] < 0))
ncut2 = np.where((newt['gb_galex'] > 0) & (newt['gb_galex'] < 5))

f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})

ax1.set_title('3", J < 13.5, VPHAS, -5 < gb < 0')
ax1.set_xlim((-0.5,2.5))
ax1.set_ylim((0,14))
ax1.set_xlabel('J - K')
ax1.set_ylabel('NUV - J')
#ax1.plot([-1.0,1.5],[-1.5,14])
ax1.axhline(y=4,color='black')

#a3 = ax1.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

a1 = ax1.scatter(star['j'][scut]-star['k'][scut],star['nuv'][scut]-star['j'][scut],edgecolor='none',alpha=0.3)
a2 = ax1.scatter(newt['j'][ncut]-newt['k'][ncut],newt['nuv_mag'][ncut]-newt['j'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)

ax1.legend([a1,a2],['Sextractor','GALEX'],scatterpoints=1,loc=2)


ax2.set_title('3", J < 13.5, VPHAS, 5 < gb < 10')
ax2.set_xlim((-0.5,2.5))
ax2.set_ylim((0,14))
ax2.set_xlabel('J - K')
ax2.set_ylabel('NUV - J')
#ax2.plot([-1.0,1.5],[-1.5,14])
ax2.axhline(y=4,color='black')
#b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

b1 = ax2.scatter(star['j'][scut2]-star['k'][scut2],star['nuv'][scut2]-star['j'][scut2],edgecolor='none',alpha=0.3)
b2 = ax2.scatter(newt['j'][ncut2]-newt['k'][ncut2],newt['nuv_mag'][ncut2]-newt['j'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)
ax2.legend([b1,b2],['Sextractor','GALEX'],scatterpoints=1,loc=2)
plt.show()



####################################################################
star = Table.read('newfield_2mass_t2_jlim_13.5_3arcsec.txt',format='ascii')
newt = Table.read('galex0data_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags.txt', format='ascii')

gbrange = 5.
scut = np.where((star['gb_sex'] > -10) & (star['gb_sex'] < -5))
scut2 = np.where((star['gb_sex'] > -5) & (star['gb_sex'] < 0))
scut3 = np.where((star['gb_sex'] > 0) & (star['gb_sex'] < 5))
scut4 = np.where((star['gb_sex'] > 5) & (star['gb_sex'] < 10))

ncut = np.where((newt['gb_galex'] > -10) & (newt['gb_galex'] < -5))
ncut2 = np.where((newt['gb_galex'] > -5) & (newt['gb_galex'] < 0))
ncut3 = np.where((newt['gb_galex'] > 0) & (newt['gb_galex'] < 5))
ncut4 = np.where((newt['gb_galex'] > 5) & (newt['gb_galex'] < 10))


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
plt.rc('legend',**{'fontsize':15})

ax1.set_title('3", J < 13.5, DF, -10 > gb > -5')
ax1.set_xlim((-1.0,2))
ax1.set_ylim((-1.5,14))
ax1.set_xlabel('J - K')
ax1.set_ylabel('NUV - J')
ax1.plot([-1.0,1.5],[-1.5,14])
a2 = ax1.scatter(newt['j'][ncut]-newt['k'][ncut],newt['nuv_mag'][ncut]-newt['j'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)
a3 = ax1.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
a1 = ax1.scatter(star['j'][scut]-star['k'][scut],star['nuv'][scut]-star['j'][scut],edgecolor='none',alpha=0.3)
ax1.legend([a1,a2,a3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)


ax2.set_title('3", J < 13.5, DF, -5 > gb > 0')
ax2.set_xlim((-1.0,2))
ax2.set_ylim((-1.5,14))
ax2.set_xlabel('J - K')
ax2.set_ylabel('NUV - J')
ax2.plot([-1.0,1.5],[-1.5,14])
b2 = ax2.scatter(newt['j'][ncut2]-newt['k'][ncut2],newt['nuv_mag'][ncut2]-newt['j'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)
b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
b1 = ax2.scatter(star['j'][scut2]-star['k'][scut2],star['nuv'][scut2]-star['j'][scut2],edgecolor='none',alpha=0.3)
ax2.legend([b1,b2,b3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)


ax3.set_title('3", J < 13.5, DF, 0 > gb > 5')
ax3.set_xlim((-1.0,2))
ax3.set_ylim((-1.5,14))
ax3.set_xlabel('J - K')
ax3.set_ylabel('NUV - J')
ax3.plot([-1.0,1.5],[-1.5,14])
c2 = ax3.scatter(newt['j'][ncut3]-newt['k'][ncut3],newt['nuv_mag'][ncut3]-newt['j'][ncut3],facecolor='none',edgecolor='red',s=30,alpha=0.3)
c3 = ax3.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax3.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
c1 = ax3.scatter(star['j'][scut3]-star['k'][scut3],star['nuv'][scut3]-star['j'][scut3],edgecolor='none',alpha=0.3)
ax3.legend([c1,c2,c3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)

ax4.set_title('3", J < 13.5, DF, 5 > gb > 10')
ax4.set_xlim((-1.0,2))
ax4.set_ylim((-1.5,14))
ax4.set_xlabel('J - K')
ax4.set_ylabel('NUV - J')
ax4.plot([-1.0,1.5],[-1.5,14])
d2 = ax4.scatter(newt['j'][ncut4]-newt['k'][ncut4],newt['nuv_mag'][ncut4]-newt['j'][ncut4],facecolor='none',edgecolor='red',s=30,alpha=0.3)
d3 = ax4.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax4.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
d1 = ax4.scatter(star['j'][scut4]-star['k'][scut4],star['nuv'][scut4]-star['j'][scut4],edgecolor='none',alpha=0.3)
ax4.legend([d1,d2,d3],['Sextractor','GALEX','Pickles'],scatterpoints=1,loc=2)


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
# NUV comparison between GALEX and SExtractor
####################################################################
star = Table.read('newfield_2mass_t2_jlim_13.5_3arcsec.txt',format='ascii')
galex = Table.read('galex0data_2mass_t2.txt',format='ascii')
pickles = Table.read('picklemags.txt',format='ascii')

stargal = SkyCoord(star.gl*u.deg,star.gb*u.deg,frame='galactic')
galexgal = SkyCoord(galex['gl_galex']*u.deg,galex['gb_galex']*u.deg,frame='galactic')
starind, galexind, angsep, dist3d = search_around_sky(stargal, galexgal, 6.*u.arcsec)
star2 = star[starind]
galex2 = galex[galexind]

gbrange = 5.
scut = np.where((np.abs(star['gb_sex']) > gbrange))
scut2 = np.where((np.abs(star['gb_sex']) < gbrange))
ncut = np.where(np.abs(newt['gb_galex']) > gbrange)
ncut2 = np.where(np.abs(newt['gb_galex']) < gbrange)
f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})
a1 = ax1.scatter(galex2['nuv_mag'][ncut],star2['nuv'][scut]-galex2['nuv_mag'][ncut],edgecolor='none',alpha=0.2)
ax1.set_title('gb > '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax1.set_xlim((13,20))
ax1.set_ylim((-2,3))
ax1.set_xlabel('GALEX NUV')
ax1.set_ylabel('SExtractor - GALEX NUV')


a2 = ax2.scatter(galex2['nuv_mag'][ncut2],star2['nuv'][scut2]-galex2['nuv_mag'][ncut2],edgecolor='none',alpha=0.2)
ax2.set_title('gb < '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax2.set_xlim((13,20))
ax2.set_ylim((-2,3))
ax2.set_xlabel('GALEX NUV')
ax2.set_ylabel('SExtractor - GALEX NUV')
plt.show()


plt.scatter(g2['nuv_mag'],s2['nuv']-1,g2['nuv_mag'],edgecolor='none',alpha=0.3)
plt.plot([12,20],[12,20],c='black')
plt.xlabel('SExtractor NUV')
plt.ylabel('GALEX NUV')
plt.title('GALEX0 vs SExtractor NUV')
plt.xlim((12,20))
plt.ylim((12,20))
plt.show()


plt.scatter(g2['nuv_mag'],s2['nuv']-g2['nuv_mag'],edgecolor='none',alpha=0.3)
plt.xlabel('GALEX NUV')
plt.ylabel('SExtractor - GALEX NUV')
plt.title('NUV Comparison')
plt.xlim((13,21))
plt.ylim((-2,2))
plt.plot([13,21],[1,1],c='black')
plt.show()

####################################################################
# Plot data w/ J limits
####################################################################
staroutcut = np.where(((star['j']-star['k']) >  0.6) & ((star['nuv']-star['j']) < 8))
zoutcut = np.where(((z['j']-z['k']) >  0.6) & ((z['nuv']-z['j']) < 8))

z = Table.read('newfield_ipac_2mass_matches_3arcsec.txt',format='ascii')
zoutcut = np.where(((z['j']-z['k']) >  0.6) & ((z['nuv']-z['j']) < 8))

img = fits.open('newfield/count_05-68_gPr_cata_10_corr.fits')[0].data
plt.imshow(img,vmin=0,vmax=1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray)
plt.scatter(z['x_new'],z['y_new'],facecolor='none',edgecolor='green',s=20,label='All')
plt.scatter(z['x_new'][zoutcut],z['y_new'][zoutcut],facecolor='none',edgecolor='red',s=20,label='J-K > 0.6')
plt.legend(loc=2,scatterpoints=1)
plt.title('Deep field, J-K > 0.6 outliers')
plt.show()


####################################################################
# Import different Jlim and arcsec search radii from 2MASS
####################################################################
from astropy.io import ascii
a2125 = Table.read('newfield_ipac_2mass_matches_2arcsec_J_lt12.5.txt',format='ascii')
a213 = Table.read('newfield_ipac_2mass_matches_2arcsec_J_lt13.txt',format='ascii')
a2135 = Table.read('newfield_ipac_2mass_matches_2arcsec_J_lt13.5.txt',format='ascii')
a214 = Table.read('newfield_ipac_2mass_matches_2arcsec_J_lt14.txt',format='ascii')
a2145 = Table.read('newfield_ipac_2mass_matches_2arcsec_J_lt14.5.txt',format='ascii')
b3125 = Table.read('newfield_ipac_2mass_matches_3arcsec_J_lt12.5.txt',format='ascii')
b313 = Table.read('newfield_ipac_2mass_matches_3arcsec_J_lt13.txt',format='ascii')
b3135 = Table.read('newfield_ipac_2mass_matches_3arcsec_J_lt13.5.txt',format='ascii')
b314 = Table.read('newfield_ipac_2mass_matches_3arcsec_J_lt14.txt',format='ascii')
b3145 = Table.read('newfield_ipac_2mass_matches_3arcsec_J_lt14.5.txt',format='ascii')
c4125 = Table.read('newfield_ipac_2mass_matches_4arcsec_J_lt12.5.txt',format='ascii')
c413 = Table.read('newfield_ipac_2mass_matches_4arcsec_J_lt13.txt',format='ascii')
c4135 = Table.read('newfield_ipac_2mass_matches_4arcsec_J_lt13.5.txt',format='ascii')
c414 = Table.read('newfield_ipac_2mass_matches_4arcsec_J_lt14.txt',format='ascii')
c4145 = Table.read('newfield_ipac_2mass_matches_4arcsec_J_lt14.5.txt',format='ascii')


#Plot radius search vs J mag search limits 
files = [a2125,a213,a2135,a214,a2145,b3125,b313,b3135,b314,b3145,c4125,c413,c4135,c414,c4145]
jlim = [12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5]
arcsec = [2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]
dx = []
for i in files:
    dx.append(np.median(i['dist_x']))

plt.scatter(jlim,dx,c=arcsec,edgecolor='none',s=40)
plt.xlabel('Jlim [mag]')
plt.ylabel('dist_x [arcsec]')
cm = plt.colorbar()
cm.set_label('Limiting radius [arcsec]')
for i in range(len(files)):
    plt.annotate(str(len(files[i])),xy=(jlim[i]+0.05,dx[i]))
plt.show()


# Field density plot vs J mag search limits
field = (6*5.5)*3*3600./2
denlo,denhi = [],[]
for i in files:
    denlo.append(len(i[np.where(np.abs(i['gb_sex']) < 5.)])/field)
    denhi.append(len(i[np.where(np.abs(i['gb_sex']) > 5.)])/field)

plt.scatter(jlim,denlo,c=arcsec,edgecolor='none',s=80,marker='s',label='b < 5')
plt.scatter(jlim,denhi,c=arcsec,edgecolor='none',s=80,marker='o',label='b > 5')
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


# NUV vs J
for i in files:
    plt.scatter(i['j'],i['nuv'],edgecolor='none',alpha=0.2)
    plt.show()

x = 0
for i in files:
    lim = np.ones(len(i))*jlim[x]
    print jlim[x]
    plt.scatter(lim,i['dist_x'],edgecolor='none',alpha=0.01)
    x+=1

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
# VPHAS matching
####################################################################
sex = Table.read('newfield_gal_ipac_tests.txt',format='ascii')
vphas = Table.read('vphas_gl_0_to_7.fits',format='fits')

vpgal = SkyCoord(vphas['RAJ2000']*u.deg,vphas['DEJ2000']*u.deg,frame='icrs')
sexind, vpind, angsep, dist3d = search_around_sky(sexgal,vpgal,5*u.arcsec)
plt.hist(angsep*3600),plt.title('SEx+VPHAS'),plt.show()
sex2 = sex[sexind]
vp2 = vphas[vpind]

sv = hstack([sex2,vp2])

sex2mass = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt13.5_tests.txt',format='ascii')
sex2massgal = SkyCoord(sex2mass['ra_2mass']*u.deg,sex2mass['dec_2mass']*u.deg,frame='icrs')
sex2mind, vp2mind, angsep, dist3d = search_around_sky(sex2massgal,vpgal,1*u.arcsec)
plt.hist(angsep*3600),plt.title('SEx+2MASS+VPHAS'),plt.show()
s2m2 = sex2mass[sex2mind]
vpm2 = vphas[vp2mind]

s2mv = hstack([s2m2,vpm2])



'''
sex = Sextractor
sv = sextractor + vphas
sex2mass = sextractor + 2MASS
s2mv = sextractor + 2mass + vphas
'''



plt.scatter(sv['g_AB']-sv['i_AB'],sv['nuv']-sv['g_AB'],edgecolor='none',alpha=0.3,label='VPHAS')
plt.scatter(vs2m['g_AB']-vs2m['i_AB'],vs2m['nuv']-vs2m['g_AB'],edgecolor='red',facecolor='none',alpha=0.3,label='VPHAS+2MASS')
plt.xlabel('g - i (ABmag)')
plt.ylabel('NUV - g (ABmag)')
plt.title('VPHAS only vs 2MASS+VPHAS ')
plt.xlim((-0.5,2))
plt.ylim((2.5,6.5))
plt.gca().invert_yaxis()
plt.legend(loc=2,scatterpoints=1)
plt.show()


plt.scatter(vgalex['g_AB']-vgalex['i_AB'],vgalex['nuv_mag']-vgalex['g_AB'],edgecolor='none',alpha=0.3,label='GALEX+VPHAS')
plt.scatter(sv['g_AB']-sv['i_AB'],sv['nuv']-sv['g_AB'],edgecolor='red',facecolor='none',alpha=0.3,label='SExtractor+VPHAS')
plt.xlabel('g - i (ABmag)')
plt.ylabel('NUV - g (ABmag)')
plt.title('VPHAS only vs 2MASS+VPHAS ')
plt.xlim((-0.5,2))
plt.ylim((2.5,6.5))
plt.gca().invert_yaxis()
plt.legend(loc=2,scatterpoints=1)
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
#plt.xlim((-0.5,8))
#plt.ylim((-10.5,11))
plt.gca().invert_xaxis()
plt.show()




plt.hist(sv['nuv'],bins=bins,label='VPHAS+SExtractor',color='green')
plt.hist(sex2mass['nuv'][s2masscut],bins=bins,label='2MASS+SExtractor',alpha=0.7,color='red')
plt.title('VPHAS vs 2MASS for abs(gb) < 4')
plt.xlabel('NUV')
plt.legend(loc=2)
plt.show()




####################################################################
# Comparing 
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

