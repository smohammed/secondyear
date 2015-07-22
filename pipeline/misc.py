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
newt = Table.read('galex120_2mass_t2.txt',format='ascii')
pickles = Table.read('picklemags.txt',format='ascii')

i = 120
gbrange = 5.
scut = np.where((np.abs(star.gb) > gbrange) & (star.gl > i) & (star.gl < i+20))
scut2 = np.where((np.abs(star.gb) < gbrange) & (star.gl > i) & (star.gl < i+20))
ncut = np.where(np.abs(newt['gb']) > gbrange)
ncut2 = np.where(np.abs(newt['gb']) < gbrange)
f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})
a1 = ax1.scatter(star.j[scut]-star.k[scut],star.nuv[scut]-star.j[scut],edgecolor='none',alpha=0.2)
ax1.set_title('gb > '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax1.set_xlim((-0.5,1.5))
ax1.set_ylim((3,12.5))
ax1.set_xlabel('J - K')
ax1.set_ylabel('NUV - J')
ax1.plot([-0.5,1.5],[3,12])
a2 = ax1.scatter(newt['j_m'][ncut]-newt['k_m'][ncut],newt['nuv_mag'][ncut]-newt['j_m'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.2)
a3 = ax1.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
ax1.legend([a1,a2,a3],['sextractor','galex','pickles'],scatterpoints=1,loc=2)


b1 = ax2.scatter(star.j[scut2]-star.k[scut2],star.nuv[scut2]-star.j[scut2],edgecolor='none',alpha=0.2)
ax2.set_title('gb < '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax2.set_xlim((-0.5,1.5))
ax2.set_ylim((3,12.5))
ax2.set_xlabel('J - K')
ax2.set_ylabel('NUV - J')
ax2.plot([-0.5,1.5],[3,12])
b2 = ax2.scatter(newt['j_m'][ncut2]-newt['k_m'][ncut2],newt['nuv_mag'][ncut2]-newt['j_m'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.2)

b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

ax2.legend([b1,b2,b3],['sextractor','galex','pickles'],scatterpoints=1,loc=2)
plt.show()
    




p1 = plt.scatter(star.j[scut]-star.k[scut],star.nuv[scut]-star.j[scut],edgecolor='none',alpha=0.3)
p2 = plt.scatter(newt['j_m']-newt['k_m'],newt['nuv']-newt['j_m'],facecolor='none',edgecolor='red',s=30,alpha=0.5)
plt.xlim((-0.5,1.5))
plt.ylim((3,12.5))
plt.xlabel('J - K')
plt.ylabel('NUV - J')
plt.title('120 < gl < 140')
plt.plot([-0.5,1.5],[3,12])
p3 = plt.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    plt.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

plt.legend([p1,p2,p3],['sextractor','GALEX','pickles'],scatterpoints=1,loc=2)
plt.show()




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
#ax1.plot([-0.5,1.5],[3,12])
#a2 = ax1.scatter(galex['j_m'][ncut]-galex['k_m'][ncut],galex['nuv_mag'][ncut]-galex['j_m'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.2)
#a3 = ax1.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
#ax1.legend([a1,a2,a3],['sextractor','galex','pickles'],scatterpoints=1,loc=2)


a2 = ax2.scatter(galex2['nuv_mag'][ncut2],star2.nuv[scut2]-galex2['nuv_mag'][ncut2],edgecolor='none',alpha=0.2)
ax2.set_title('gb < '+str(gbrange)+', gl = '+str(i)+' to '+str(i+20))
ax2.set_xlim((13,20))
ax2.set_ylim((-2,3))
ax2.set_xlabel('GALEX NUV')
ax2.set_ylabel('SExtractor - GALEX NUV')
#ax2.plot([-0.5,1.5],[3,12])
#b2 = ax2.scatter(galex['j_m'][ncut2]-galex['k_m'][ncut2],galex['nuv_mag'][ncut2]-galex['j_m'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.2)

#b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
#for j in range(len(pickles)):
#    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

#ax2.legend([b1,b2,b3],['sextractor','galex','pickles'],scatterpoints=1,loc=2)
plt.show()
    