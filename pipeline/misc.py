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


pickles = Table.read('picklemags.txt',format='ascii')

i = 120
scut = np.where((np.abs(star.gb > 5.)) & (star.gl > i) & (star.gl < i+20))
scut2 = np.where((np.abs(star.gb < 5.)) & (star.gl > i) & (star.gl < i+20))
ncut = np.where(np.abs(newt['gb'] > 5))
ncut2 = np.where(np.abs(newt['gb'] < 5))
f, (ax1, ax2) = plt.subplots(1, 2)
plt.rc('legend',**{'fontsize':15})
a1 = ax1.scatter(star.j[scut]-star.k[scut],star.nuv[scut]-star.j[scut],edgecolor='none',alpha=0.2)
ax1.set_title('gb > 5, gl = '+str(i)+' to '+str(i+20))
ax1.set_xlim((-0.5,1.5))
ax1.set_ylim((3,12.5))
ax1.set_xlabel('J - K')
ax1.set_ylabel('NUV - J')
ax1.plot([-0.5,1.5],[3,12])
a2 = ax1.scatter(newt['j_m'][ncut]-newt['k_m'][ncut],newt['nuv'][ncut]-newt['j_m'][ncut],facecolor='none',edgecolor='red',s=30,alpha=0.5)
a3 = ax1.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax1.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)
ax1.legend([a1,a2,a3],['sextractor','galex','pickles'],scatterpoints=1,loc=2)


b1 = ax2.scatter(star.j[scut2]-star.k[scut2],star.nuv[scut2]-star.j[scut2],edgecolor='none',alpha=0.2)
ax2.set_title('gb < 5, gl = '+str(i)+' to '+str(i+20))
ax2.set_xlim((-0.5,1.5))
ax2.set_ylim((3,12.5))
ax2.set_xlabel('J - K')
ax2.set_ylabel('NUV - J')
ax2.plot([-0.5,1.5],[3,12])
b2 = ax2.scatter(newt['j_m'][ncut2]-newt['k_m'][ncut2],newt['nuv'][ncut2]-newt['j_m'][ncut2],facecolor='none',edgecolor='red',s=30,alpha=0.5)

b3 = ax2.scatter(pickles['j']-pickles['k'],pickles['nuv']-pickles['j'],c='black',s=5)
for j in range(len(pickles)):
    ax2.annotate(pickles['name'][j][:-4],xy=(pickles['j'][j]-pickles['k'][j],pickles['nuv'][j]-pickles['j'][j]),size=12)

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





