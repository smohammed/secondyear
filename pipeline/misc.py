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

##########################################################
# g vs g-r
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

# B - V version
plt.scatter(c2['B_AB']-c2['V_AB'], c2['FE_H'], edgecolor='none', c=c2['ALPHAFE'], s=80, vmin=-0.05, vmax=0.3, **{"zorder":100})
plt.errorbar(c2['B_AB']-c2['V_AB'], c2['FE_H'], xerr=c2['B_ABerr']-c2['V_ABerr'], yerr=c2['FE_H_ERR'], fmt=None, marker=None, mew=0, **{"zorder":0})

#plt.xlim((-1, 0.4))
plt.xlabel('B - V')
#plt.xlabel('(B - E$_{B-V}$ * 3.626) - (V - E$_{B-V}$ * 2.742)')
plt.ylabel('Fe/H')
#plt.title('GAIS + Bovy RC stars, no extinction')
plt.colorbar().set_label('Alpha/Fe')
plt.show()

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


#################################################
# Combine all starcat files
#################################################
scans = ['18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '189-199', '198-208', '207-217', '216-226', '270-280', '279-289', '288-298', '297-307', '351-1', '72-82', '144-154', '225-235', '306-316', '81-91', '153-163', '234-244', '315-325', '90-100', '162-172', '243-253', '324-334', '9-19', '171-181', '252-262', '333-343', '99-109', '180-190', '261-271', '342-352']

alldata = Table.read('starcat_0-10_03_26_2018.txt', format='ascii')

for region in scans:
	a = Table.read('starcat_'+region+'_03_26_2018.txt', format='ascii')
	a = a[np.where((a['FWHM_IMAGE'] > 0) & (a['FWHM_IMAGE'] < 10))]
	alldata = vstack([alldata, a])
	print region

ascii.write(alldata, 'starcat_allscans_03-26-18.txt', format='basic', overwrite=True)


pscats = ['0-100', '100-200', '200-300', '300-360']

for scan in pscats:
    ps = fits.open('ps1/ps1_planecut_'+scan+'_g10-20.fits')[1].data
    psgal = SkyCoord(ps['glmean']*u.deg, ps['gbmean']*u.deg, frame='galactic')
    catgal = SkyCoord(alldata['gl']*u.deg,alldata['gb']*u.deg, frame='galactic')
    catind, psind, angsep, ang3d = search_around_sky(catgal, psgal, 3*u.arcsec)
    comb = hstack([Table(alldata[catind]), Table(ps[psind])])
    comb['angsep'] = angsep
    ascii.write('starcat_ps1'+str(scan)+'_g10-20.txt', format='basic')


######################################################################
# Color vs TEFF vs Fe/H
######################################################################
m = (0.095-.21)/(0+0.8)
b = 0.095
thick, = np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin, = np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))
#y1 = (rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742)
y1 = rc['phot_bp_mean_mag']-rc['phot_rp_mean_mag']
y2 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
x = rc['TEFF']
m1, b1, rval1, pval1, stderr1 = stats.linregress(x, y1)
line1 = m1*x + b1
err1 = np.sqrt(np.sum((line1-y1)**2/len(y1)))
m2, b2, rval2, pval2, stderr2 = stats.linregress(x, y2)
line2 = m2*x + b2
err2 = np.sqrt(np.sum((line2-y2)**2/len(y2)))

fig, (ax1, ax2) = plt.subplots(2, 1)
cmap = ax1.scatter(x[thin], y1[thin], s=80, c=rc['FE_H'][thin], cmap=cm.viridis_r, vmin=-.5, vmax=.35, marker='D', label='Thin disk', edgecolor='black')
ax1.scatter(x[thick], y1[thick], s=80, c=rc['FE_H'][thick], cmap=cm.viridis_r, vmin=-.5, vmax=.35, marker='s', label='Thick disk', edgecolor='black', linewidth=3)
#ax1.errorbar(x[thin], y1[thin], xerr=rc['TEFF_ERR'][thin], yerr=(rc['Berr_apass']-rc['Verr_apass'])[thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax1.errorbar(x[thin], y1[thin], xerr=rc['TEFF_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax1.plot(x, line1, linewidth=5, c='red', zorder=10)

ax2.scatter(x[thin], y2[thin], s=80, c=rc['FE_H'][thin], cmap=cm.viridis_r, vmin=-.5, vmax=.35, marker='D', edgecolor='black')
ax2.scatter(x[thick], y2[thick], s=80, c=rc['FE_H'][thick], cmap=cm.viridis_r, vmin=-.5, vmax=.35, marker='s', edgecolor='black', linewidth=3)
#ax2.errorbar(x[thin], y2[thin], xerr=rc['TEFF_ERR'][thin], yerr=(rc['nuv_magerr']-rc['Gerr'])[thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax2.errorbar(x[thin], y2[thin], xerr=rc['TEFF_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax2.plot(x, line2, linewidth=5, c='red', zorder=10)

#ax1.set_xlabel('T$_{eff}$')
ax1.set_ylabel('BP - RP')
ax2.set_xlabel('T$_{eff}$')
ax2.set_ylabel('(NUV - G)$_0$')
#ax1.set_xlim((4300, 5200))
#ax2.set_xlim((4300, 5200))
#ax1.set_ylim((0.5, 1.6))
#ax2.set_ylim((5, 12))
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


##########################################
# Alpha/Fe vs Fe/H vs NUV - G
##########################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0))
rc = rc[q]

m = (0.095-.21)/(0+0.8)
b = 0.095
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))
nuvg = ((rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303))

#afe_apo = np.where(rc['ALPHA_M'] > 0)
#rc['ALPHAFE'][afe_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FE_H'])[afe_apo]
#alphafe = rc['ALPHAFE']

# for now use this. err is only for apo
#afeerr = (alphafe * np.sqrt((rc['ALPHA_M_ERR']/rc['ALPHA_M'])**2 + (rc['M_H_ERR']/rc['M_H'])**2 + (rc['FE_H_ERR']/rc['FE_H'])**2))

plt.scatter(rc['FE_H'][thin], rc['ALPHAFE'][thin], s=20, edgecolor='none', c=nuvg[thin], label='Thin disk', vmin=7, vmax=11, marker='D', cmap=cm.jet)

plt.scatter(rc['FE_H'][thick], rc['ALPHAFE'][thick], c=nuvg[thick], s=20, marker='s', label='Thick disk', vmin=7, vmax=11, linewidth=1, cmap=cm.jet)

#plt.errorbar(rc['FE_H'], rc['ALPHAFE'], xerr=(rc['FE_H_err']), yerr=rc['ALPHAFE_ERR'], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
plt.errorbar(rc['Fe_H'], rc['ALPHAFE'], xerr=(rc['FE_H_err']), ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})


#plt.scatter(rc['FE_H'][thin], rc['ALPHAFE'][thin], s=20, edgecolor='none', label='Thin disk', vmin=7, vmax=11, marker='D')

#plt.scatter(rc['FE_H'][thick], rc['ALPHAFE'][thick], s=20, marker='s', label='Thick disk', vmin=7, vmax=11, linewidth=1)


plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha$/Fe]')
#plt.xlim((-0.7, 0.5))
#plt.ylim((-0.05, 0.3))
plt.colorbar().set_label('(NUV - G)$_0$')
leg = plt.legend(scatterpoints=1)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
plt.show()

#################################################
# rc fit lines for Fe/H vs NUV - G vs Alpha/Fe
#################################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0))
rc = rc[q]

m = (0.095-.21)/(0+0.8)
b = 0.095
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

# Key:
# a is not dust corrected. b is.
# 1 includes the outliers/possible binaries. 2 removes them

xa1 = rc['nuv_mag']-rc['phot_g_mean_mag']
xb1 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*3.303)
y1 = rc['FE_H']

q = np.where(xb1 > 5)
w = np.where((xb1 < 5) & (rc['Fe_H'] > -0.6))
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


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

alpha = 1
# Plotting with no dust correction first
cmap = ax1.scatter(xa1[thick], y1[thick], c=rc['ALPHAFE'][thick], s=80, vmin=-0.05, vmax=0.3, marker='s', edgecolor='blue', linewidth=1, alpha=alpha, cmap=cm.plasma_r, **{"zorder":5})
#ax1.scatter(xa1[thin], y1[thin], c=rc['ALPHAFE'][thin], s=80, vmin=-0.05, vmax=0.3, marker='D', alpha=alpha, edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
#ax1.errorbar(xa1[thin], y1[thin], xerr=(rc['nuv_magerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax1.errorbar(xa1[thick], y1[thick], xerr=(rc['nuv_magerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax1.plot(xp, pa1(xp), linewidth=4, c='red', zorder=10) # no dust, all points
ax1.plot(xp, pa2(xp), linewidth=4, c='black', zorder=10) # no dust, only black points
ax1.scatter(xa1[w], y1[w], c=rc['ALPHAFE'][w], s=80, vmin=-0.05, vmax=0.3, marker='D', cmap=cm.plasma_r, edgecolor='red', linewidth=2, **{"zorder":5})


ax2.scatter(xb1[thick], y1[thick], c=rc['ALPHAFE'][thick], s=80, label='Thick disk', vmin=-0.05, vmax=0.3, marker='s', edgecolor='blue', linewidth=1, cmap=cm.plasma_r, **{"zorder":5})
#ax2.scatter(xb1[thin], y1[thin], c=rc['ALPHAFE'][thin], s=80, label='Thin disk', vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
#ax2.errorbar(xb1[thin], y1[thin], xerr=(rc['nuv_magerr'])[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})

ax2.errorbar(xb1[thick], y1[thick], xerr=(rc['nuv_magerr'])[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt=None, marker=None, mew=0, elinewidth=1.3, **{"zorder":0})
ax2.plot(xp, pb1(xp), linewidth=4, c='red', zorder=10)
ax2.plot(xp, pb2(xp), linewidth=4, c='black', zorder=10)
ax2.scatter(xb1[w], y1[w], c=rc['ALPHAFE'][w], s=80, vmin=-0.05, vmax=0.3, marker='D', cmap=cm.plasma_r, edgecolor='red', linewidth=2, **{"zorder":5})

ax1.set_xlabel('NUV - G')
ax2.set_xlabel('(NUV - G)$_{0}$')
ax1.set_ylabel('[Fe/H]')
leg = ax2.legend(scatterpoints=1, loc=4)
leg.legendHandles[0].set_color('black')
#leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
#leg.legendHandles[1]._sizes = [150]
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
ax1.annotate('[Fe/H] = 0.25 * (NUV-G) - 2.29', xy=(2, 0.69), color='red', size=15)
ax1.annotate('[Fe/H] = 0.25 * (NUV-G) - 2.26 ', xy=(2, 0.6), color='black', size=15)
ax2.annotate('[Fe/H] = 0.26 * (NUV-G)$_0$ - 2.29', xy=(1.3, 0.69), color='red', size=15)
ax2.annotate('[Fe/H] = 0.27 * (NUV-G)$_0$ - 2.33', xy=(1.3, 0.6), color='black', size=15)

ax1.annotate('$\sigma$ = 0.018', xy=(2, 0.53), color='red', size=15) # a1
ax1.annotate('$\sigma$ = 0.018', xy=(2, 0.45), color='black', size=15) # a2
ax2.annotate('$\sigma$ = 0.018', xy=(1.3, 0.53), color='red', size=15) # b1
ax2.annotate('$\sigma$ = 0.018', xy=(1.3, 0.45), color='black', size=15) # b2
plt.suptitle('Only thick disk')
fig.colorbar(cmap, cax=cbar_ax).set_label(r'[$\alpha$/Fe]')
plt.show()



######################################################################
# Add cmd lines to cmd
######################################################################
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

######################################################################
# Combine and format RC table
######################################################################
sg = fits.open('sex_tgas_dust_interp_02-28-2018.fits')[1].data
sg = sg[~np.isnan(sg['ebv'])]
sggal = SkyCoord(sg['gl_sex']*u.deg, sg['gb_sex']*u.deg, frame='galactic')

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
bv = (rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742)
Berr = rc['Berr_apass']
bverr = np.sqrt(rc['Berr_apass']**2 - rc['Verr_apass']**2)
bverr[np.isnan(bverr)] = Berr[np.isnan(bverr)]
ebv = rc['ebv']
dm = rc['distmod']
feh = rc['FE_H']
feherr = rc['FE_H_ERR']
teff = rc['TEFF']
alphafe = rc['ALPHAFE']
afeerr = rc['ALPHAFE_ERR']

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
	hdu = fits.open('count_map_'+scan+'_in.fits')[0]
	wcs = WCS(hdu.header)
	bkgd = fits.open('background_im1_'+scan+'.fits')[0].data
	t = Table.read('starcat_'+scan+'_03_26_2018.txt', format='ascii')
	t = t[np.where((t['FWHM_IMAGE'] > 0) & (t['FWHM_IMAGE'] < 10) & (t['expsum'] > 10.))]

	if scan == '0-10':
		t['gl'][np.where(t['gl'] > 355)] = t['gl'][np.where(t['gl'] > 355)] - 360.

	print(scan)

	fig = plt.figure()
	fig.add_subplot(221, projection=wcs)
	plt.imshow(hdu.data, origin='lower', cmap=cm.gray, aspect='auto', vmin=0, vmax=0.1)
	plt.title(scan)
	
	fig.add_subplot(222, projection=wcs)
	plt.imshow(bkgd, origin='lower', cmap=cm.gray, aspect='auto', vmin=0, vmax=0.1)

	fig.add_subplot(223)
	cbar = plt.scatter(t['gl'], t['gb'], c=t['FWHM_IMAGE'], vmin=0, vmax=10, s=3)
	plt.colorbar(cbar).set_label('FWHM')
	plt.gca().invert_xaxis()
	plt.xlabel('gl')
	plt.ylabel('gb')
	plt.title(scan+', 0 > FWHM > 10')

	fig.add_subplot(224)
	plt.scatter(t['nuv'], t['FWHM_IMAGE'], s=1, c='black')
	plt.ylim(0, 10)
	plt.title('N = '+str(len(t)))

	plt.savefig('03-26-2018_'+scan+'_4panel_scanplots.png')
	plt.clf()
	del bkgd
	del hdu

#############################################################
# Get values from maps
#############################################################
scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '189-199', '198-208', '207-217', '216-226', '270-280', '279-289', '288-298', '297-307', '351-1', '72-82', '144-154', '225-235', '306-316', '81-91', '153-163', '234-244', '315-325', '90-100', '162-172', '243-253', '324-334', '9-19', '171-181', '252-262', '333-343', '99-109', '180-190', '261-271', '342-352']

from photutils import aperture_photometry
from photutils import CircularAperture

for scan in scans:
	print(scan)
	cat = Table.read('starcat_'+scan+'_03_26_2018.txt', format='ascii')
	bkgd = fits.open('background_im1_'+scan+'.fits')[0].data
	positions = [cat['X_IMAGE'], cat['Y_IMAGE']]
	apertures = CircularAperture(positions, r=3.)
	sums = aperture_photometry(bkgd, apertures)
	cat['bkgdsum'] = sums['aperture_sum']
	ascii.write(cat, 'starcat_'+scan+'_03_26_2018.txt', format='basic')
