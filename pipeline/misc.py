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

star = fits.open('starcatalog_120_rand_deadtime.fits')[1].data



#NUV - B, B - V

star = Table.read('starcatalog_05-68_2mass_t2.txt', format='ascii')
newt = Table.read('galex120_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags.txt', format='ascii')
gbrange = 5.
scut = np.where((np.abs(star['gb_sex']) > gbrange))
scut2 = np.where((np.abs(star['gb_sex']) < gbrange))
ncut = np.where(np.abs(newt['gb']) > gbrange)
ncut2 = np.where(np.abs(newt['gb']) < gbrange)

f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})

ax1.set_title('New Field, gb > '+str(gbrange))
ax1.set_xlim((0,2))
ax1.set_ylim((0,10))
ax1.set_xlabel('B - V')
ax1.set_ylabel('NUV - B')
ax1.plot([0,2],[0,10])
a2 = ax1.scatter(newt['BJmag'][ncut]-newt['VJmag'][ncut],newt['nuv_mag'][ncut]-newt['BJmag'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.3)
a3 = ax1.scatter(pickles['b']-pickles['v'],pickles['nuv']-pickles['b'],c='black',s=5)
for j in range(len(pickles)):
    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['b'][j]-pickles['v'][j],pickles['nuv'][j]-pickles['b'][j]),size=12)

a1 = ax1.scatter(star['BJmag'][scut]-star['VJmag'][scut],star['nuv'][scut]-star['BJmag'][scut],edgecolor='none',alpha=0.5)
ax1.legend([a1,a2,a3],['Sextractor','galex','pickles'],scatterpoints=1,loc=2)


ax2.set_title('New Field, gb < '+str(gbrange))
ax2.set_xlim((0,2))
ax2.set_ylim((0,10))
ax2.set_xlabel('B - V')
ax2.set_ylabel('NUV - B')
ax2.plot([0,2],[0,10])
b2 = ax2.scatter(newt['BJmag'][ncut2]-newt['VJmag'][ncut2],newt['nuv_mag'][ncut2]-newt['BJmag'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.3)
b3 = ax2.scatter(pickles['b']-pickles['v'],pickles['nuv']-pickles['b'],c='black',s=5)
for j in range(len(pickles)):
    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['b'][j]-pickles['v'][j],pickles['nuv'][j]-pickles['b'][j]),size=12)

b1 = ax2.scatter(star['BJmag'][scut2]-star['VJmag'][scut2],star['nuv'][scut2]-star['BJmag'][scut2],edgecolor='none',alpha=0.5)
ax2.legend([b1,b2,b3],['Sextractor','galex','pickles'],scatterpoints=1,loc=2)
plt.show()
    


# NUV - J vs J - K

star = Table.read('starcatalog_05-68_2mass_t2.txt', format='ascii')
newt = Table.read('galex120_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags.txt', format='ascii')

gbrange = 5.
scut = np.where((np.abs(star['gb_sex']) > gbrange))
scut2 = np.where((np.abs(star['gb_sex']) < gbrange))
ncut = np.where(np.abs(newt['gb']) > gbrange)
ncut2 = np.where(np.abs(newt['gb']) < gbrange)

f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})

ax1.set_title('3", J < 13.5, Deeper Field, gb > '+str(gbrange))
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


ax1.legend([a1,a2,a3],['Sextractor','Galex','pickles'],scatterpoints=1,loc=2)


ax2.set_title('3", J < 13.5, Deeper Field, gb < '+str(gbrange))
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
ax2.legend([b1,b2,b3],['Sextractor','Galex','pickles'],scatterpoints=1,loc=2)
plt.show()



c.rename_column('cntr_01','cntr')
c.rename_column('number_01','number')
c.rename_column('x_image_01','x_image')
c.rename_column('y_image_01','y_image')
c.rename_column('flux_auto_01','flux_auto')
c.rename_column('fluxerr_auto_01','fluxerr_auto')
c.rename_column('x_new_01','x_new')
c.rename_column('y_new_01','y_new')
c.rename_column('nuv_01','nuv')
c.rename_column('gl_01','gl_sex')
c.rename_column('gb_01','gb_sex')
c.rename_column('ra_01','ra_sex')
c.rename_column('dec_01','dec_sex')
c.rename_column('ra','ra_2mass')
c.rename_column('dec','dec_2mass')
c.rename_column('j_m','j')
c.rename_column('h_m','h')
c.rename_column('k_m','k')



weight1 = np.ones(np.shape(art)) * 100
weight2 = np.ones(np.shape(art)) * 100
for i in range(2000):
    for j in range(2000):
        try:
            if ((art[i-1][j-1]+art[i-1][j]+art[i-1][j+1]+art[i][j-1]+art[i][j+1]+art[i+1][j-1]+art[i+1][j]+art[i+1][j+1]+art[i][j]) == 0.):
                weight1[i][j] = 0.
            if ((art[i-1][j-1]+art[i-1][j]+art[i-1][j+1]+art[i][j-1]+art[i][j+1]+art[i+1][j-1]+art[i+1][j]+art[i+1][j+1]) == 0.) and (art[i][j] > 0.08):
                #weight1[i][j] = 0.
                weight1[i-1][j-1] = 100.
                weight1[i-1][j] = 100.
                weight1[i-1][j+1] = 100.
                weight1[i][j-1] = 100.
                weight1[i][j+1] = 100.
                weight1[i+1][j-1] = 100.
                weight1[i+1][j] = 100.
                weight1[i+1][j+1] = 100.
        except IndexError:
            pass


star = fits.open('starcatalog_120_rand_deadtime.fits')[1].data
galex = Table.read('galex120_2mass_t2.txt',format='ascii')
pickles = Table.read('picklemags.txt',format='ascii')

stargal = SkyCoord(star.gl*u.deg,star.gb*u.deg,frame='galactic')
galexgal = SkyCoord(galex['gl']*u.deg,galex['gb']*u.deg,frame='galactic')
starind, galexind, angsep, dist3d = search_around_sky(stargal, galexgal, 6.*u.arcsec)
star2 = star[starind]
galex2 = galex[galexind]

i = 120
gbrange = 5.
scut = np.where((np.abs(star2.gb) > gbrange) & (star2.gl > i) & (star2.gl < i+20))
scut2 = np.where((np.abs(star2.gb) < gbrange) & (star2.gl > i) & (star2.gl < i+20))
ncut = np.where(np.abs(galex2['gb']) > gbrange)
ncut2 = np.where(np.abs(galex2['gb']) < gbrange)
f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})
a1 = ax1.scatter(galex2['nuv_mag'][ncut],star2.nuv[scut]-galex2['nuv_mag'][ncut],edgecolor='none',alpha=0.2)
ax1.set_title('gb > '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax1.set_xlim((13,20))
ax1.set_ylim((-2,3))
ax1.set_xlabel('GALEX NUV')
ax1.set_ylabel('SExtractor - GALEX NUV')


a2 = ax2.scatter(galex2['nuv_mag'][ncut2],star2.nuv[scut2]-galex2['nuv_mag'][ncut2],edgecolor='none',alpha=0.2)
ax2.set_title('gb < '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax2.set_xlim((13,20))
ax2.set_ylim((-2,3))
ax2.set_xlabel('GALEX NUV')
ax2.set_ylabel('SExtractor - GALEX NUV')
plt.show()




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


files = [a2125,a213,a2135,a214,a2145,b3125,b313,b3135,b314,b3145,c4125,c413,c4135,c414,c4145]
arcsec = ['2','2','2','2','2','3','3','3','3','3','4','4','4','4','4']
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


field = (6*5.5)*3*3600.
density = []
for i in files:
    density.append(len(i)/field)


plt.scatter(jlim,density,c=arcsec,edgecolor='none',s=40)
plt.xlabel('Jlim [mag]')
plt.ylabel('Counts/area [#/arcsec$^2$]')
cm = plt.colorbar()
cm.set_label('Limiting radius [arcsec]')
for i in range(len(files)):
    plt.annotate(str(len(files[i])),xy=(jlim[i]+0.05,density[i]))
plt.title('2MASS match counts by J_lim and radius search')
plt.show()



for i in files:
    plt.scatter(i['j'],i['nuv'],edgecolor='none',alpha=0.2)
    plt.show()


x = 0
for i in files:
    lim = np.ones(len(i))*jlim[x]
    print jlim[x]
    plt.scatter(lim,i['dist_x'],edgecolor='none',alpha=0.01)
    x+=1



files = [a2125,a213,a2135,a214,a2145,b3125,b313,b3135,b314,b3145,c4125,c413,c4135,c414,c4145]

arcsec = ['2','2','2','2','2','3','3','3','3','3','4','4','4','4','4']
jlim = [12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5,12.5,13,13.5,14,14.5]
x = 0
for i in files:
    star = i
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
    ascii.write(star,'newfield_ipac_2mass_matches_'+arcsec[x]+'arcsec_J_lt'+str(jlim[x])+'.txt',format='ipac')  
    print jlim[x]
    x=x+1
    print x