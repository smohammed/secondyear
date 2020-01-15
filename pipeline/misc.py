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
# CMD with MG and MNUV
######################################################################
sg = fits.open('gais_tgas_apass_dust.fits')[1].data
sg = sg[~np.isnan(sg['ebv'])]
artefact
fig, (ax1, ax2) = plt.subplots(2, 1)artefact
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
# Combine all starcat files with no overlap
#################################################
scans = ['9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']

alldata = Table.read('starcat_0-10_06_06_2019.txt', format='ascii')
alldata = alldata[np.where(alldata['gl'] < 300.)]

for region in scans:
    a = Table.read('starcat_'+region+'_06_06_2019.txt', format='ascii')
    if region == '351-1':
        a = a[np.where(a['gl'] > 300.)]
    a = a[np.where((a['FWHM_IMAGE'] > 0) & (a['FWHM_IMAGE'] < 10) & (a['gl'] > np.min(a['gl']+1.562)))]
    alldata = vstack([alldata, a])
    print region

ascii.write(alldata, 'starcat_allscans_06-06-19.txt', format='basic', overwrite=True)

#################################################
# Match PS1 to scan data
#################################################
cat = fits.open('starcat_allscans_10-12-18.fits')[1].data
catgal = SkyCoord(cat['gl']*u.deg,cat['gb']*u.deg, frame='galactic')
pscans = ['0-100', '100-200', '200-300', '300-360']

for scan in pscans:
    print scan
    ps = fits.open('ps1/ps1_planecut_'+scan+'_g10-20.fits')[1].data
    psgal = SkyCoord(ps['glmean']*u.deg, ps['gbmean']*u.deg, frame='galactic')
    catind, psind, angsep, ang3d = search_around_sky(catgal, psgal, 3*u.arcsec)
    comb = hstack([Table(cat[catind]), Table(ps[psind])])
    comb['angsep'] = angsep
    ascii.write(comb, 'starcat_ps1'+scan+'_g10-20.txt', format='basic')
# Then use topcat to combine

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


################################################################
# Image plots for each scan output
################################################################
scans = ['0-10', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']


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
# Get values from count, exp, and bkgd map
#############################################################
scans = ['0-10', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']


for scan in scans:
    print(scan)
    cat = Table.read('starcat_'+scan+'_06_06_2019.txt', format='ascii')

    #bkgd = fits.open('background_im1_'+scan+'.fits')[0].data
    exp = fits.open('count_map_'+scan+'_exp.fits')[0].data
    #ct = fits.open('count_map_'+scan+'_count.fits')[0].data

    positions = [cat['X_IMAGE'], cat['Y_IMAGE']]
    #bkgdval = []
    expval = []
    #ctval = []
    for line in range(len(cat)):
        try:
            #ctval.append(ct[int(cat['Y_IMAGE'][line]), int(cat['X_IMAGE'][line])])
            #bkgdval.append(bkgd[int(cat['Y_IMAGE'][line]), int(cat['X_IMAGE'][line])])
            expval.append(exp[int(cat['Y_IMAGE'][line]), int(cat['X_IMAGE'][line])])
        except IndexError:
            #ctval.append(ct[int(cat['Y_IMAGE'][line])-1, int(cat['X_IMAGE'][line])-1])
            #bkgdval.append(bkgd[int(cat['Y_IMAGE'][line])-1, int(cat['X_IMAGE'][line])-1])
            expval.append(exp[int(cat['Y_IMAGE'][line])-1, int(cat['X_IMAGE'][line])-1])
    #cat['ctsum'] = ctval
    #cat['bkgdsum'] = bkgdval
    cat['expsum'] = expval

    ascii.write(cat, 'starcat_'+scan+'_06_06_2019.txt', format='basic', overwrite=True)



#############################################################
# CMD for plane survey vs GAIS
#############################################################
cg = fits.open('plane_gaiadr2_dust_09_27_18.fits')[1].data
gag = fits.open('gais_gaiadr2_negparcuts_dust_09_27_18.fits')[1].data

negpar = np.where((cg['dist'] > 0) & (cg['visibility_periods_used'] > 8) & (cg['phot_bp_mean_mag'] > 0) & (cg['phot_rp_mean_mag'] > 0) & (cg['expsum'] > 10) & (cg['ebv'] > 0) & (cg['parallax_error']/cg['parallax'] < 0.1))
#cut = np.where((gag['parallax'] > 0) & (gag['phot_g_mean_mag'] > 0) & (gag['phot_bp_mean_mag'] > 0) & (gag['phot_rp_mean_mag'] > 0) & (gag['visibility_periods_used'] > 8) & (gag['mag_nuv'] > 0) & (gag['glat'] > -10) & (gag['glat'] < 10))

cg = cg[negpar]
gag = gag[np.where((gag['ebv'] > 0) & (gag['parallax_error']/gag['parallax'] < 0.1))]

threshold = 1000
bins = 40

fig, axes = plt.subplots(nrows=2, ncols=2)

scatter_contour(gag['mag_nuv']-gag['phot_g_mean_mag'], gag['phot_g_mean_mag']-gag['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 0])

scatter_contour(cg['nuv']-cg['phot_g_mean_mag'], cg['phot_g_mean_mag']-cg['distmod'], threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-1,12], [-2,14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[0, 1])

# If I want a RC paper CMD template
#scatter_contour(cg['phot_bp_mean_mag']-cg['phot_rp_mean_mag'], cg['phot_g_mean_mag']-cg['distmod']-cg['ebv']*2.85, threshold=threshold, log_counts=True, histogram2d_args=dict(bins=bins, range=[[-0.5, 2.5],[-2, 14]]), plot_args=dict(color='k', markersize=1, alpha=0.1), contour_args=dict(cmap=cm.gray, zorder=10), ax=axes[1, 0])

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

#############################################################
# dx dy histogram by magnitude
#############################################################
cg = fits.open('plane_gaiadr2_comb_05_10_18.fits')[1].data

dx = (cg['l']-cg['gl'])*3600
dy = (cg['b']-cg['gb'])*3600

for mag in range(12, 21,2):
    cut = np.where((cg['nuv'] > mag) & (cg['nuv'] < mag+2))
    dxmag = dx[cut]
    dymag = dy[cut]
    plt.hist(np.sqrt(dxmag**2+dymag**2), histtype='step', fill=False, stacked=True, label=str(mag)+'-'+str(mag+2), range=[0,3], bins=10)

plt.legend(scatterpoints=1, loc=1)
plt.xlabel('Radius (Gaia DR2 - plane)')
plt.show()

#############################################################
# Sky plot matching the overlap with GAIS
#############################################################
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
gal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
gl = gal.galactic.l.degree
gb = gal.galactic.b.degree

gl1 = gl[np.where((gl > 0.) & (gl < 90.))]
gl2 = gl[np.where((gl > 90.) & (gl < 180.))]
gl3 = gl[np.where((gl > 180.) & (gl < 270.))]
gl4 = gl[np.where((gl > 270.) & (gl < 360.))]

gb1 = gb[np.where((gl > 0.) & (gl < 90.))]
gb2 = gb[np.where((gl > 90.) & (gl < 180.))]
gb3 = gb[np.where((gl > 180.) & (gl < 270.))]
gb4 = gb[np.where((gl > 270.) & (gl < 360.))]

cat = fits.open('starcat_allscans_03-26-18.fits')[1].data
cat1 = cat[np.where((cat['gl'] > 0.) & (cat['gl'] < 90.))]
cat2 = cat[np.where((cat['gl'] > 90.) & (cat['gl'] < 180.))]
cat3 = cat[np.where((cat['gl'] > 180.) & (cat['gl'] < 270.))]
cat4 = cat[np.where((cat['gl'] > 270.) & (cat['gl'] < 360.))]

fig, axes = plt.subplots(4, 1, sharey=True)

cbar = axes[0].scatter(cat1['gl'], cat1['gb'], c=cat1['nuv'], vmin=12, vmax=20, s=3, alpha=0.5)
axes[1].scatter(cat2['gl'], cat2['gb'], c=cat2['nuv'], vmin=12, vmax=20, s=3, alpha=0.5)
axes[2].scatter(cat3['gl'], cat3['gb'], c=cat3['nuv'], vmin=12, vmax=20, s=3, alpha=0.5)
axes[3].scatter(cat4['gl'], cat4['gb'], c=cat4['nuv'], vmin=12, vmax=20, s=3, alpha=0.5)

axes[0].scatter(gl1, gb1, facecolor='none', edgecolor='red', s=75)
axes[1].scatter(gl2, gb2, facecolor='none', edgecolor='red', s=75)
axes[2].scatter(gl3, gb3, facecolor='none', edgecolor='red', s=75)
axes[3].scatter(gl4, gb4, facecolor='none', edgecolor='red', s=75)

fig.subplots_adjust(right=0.9)
fig.colorbar(cbar, cax=fig.add_axes([0.92, 0.15, 0.02, 0.7])).set_label('NUV')

axes[0].set_xlim(0, 90)
axes[1].set_xlim(90, 180)
axes[2].set_xlim(180, 270)
axes[3].set_xlim(270, 360)

axes[3].set_xlabel('gl')
axes[0].set_ylabel('gb')
axes[1].set_ylabel('gb')
axes[2].set_ylabel('gb')
axes[3].set_ylabel('gb')
plt.show()

# Same as above but normalizing the histogram
fig, axes = plt.subplots(4, 1, sharey=True)
cbar = axes[0].hist2d(cat1['gl'], cat1['gb'], bins=(90, 20), norm=LogNorm())
axes[1].hist2d(cat2['gl'], cat2['gb'], bins=(90, 20), norm=LogNorm())
axes[2].hist2d(cat3['gl'], cat3['gb'], bins=(90, 20), norm=LogNorm())
axes[3].hist2d(cat4['gl'], cat4['gb'], bins=(90, 20), norm=LogNorm())

axes[0].scatter(gl1, gb1, facecolor='none', edgecolor='red', s=75)
axes[1].scatter(gl2, gb2, facecolor='none', edgecolor='red', s=75)
axes[2].scatter(gl3, gb3, facecolor='none', edgecolor='red', s=75)
axes[3].scatter(gl4, gb4, facecolor='none', edgecolor='red', s=75)

fig.subplots_adjust(right=0.9)
fig.colorbar(cbar[3], cax=fig.add_axes([0.92, 0.15, 0.02, 0.7]))

axes[0].set_xlim(0, 90)
axes[1].set_xlim(90, 180)
axes[2].set_xlim(180, 270)
axes[3].set_xlim(270, 360)

axes[3].set_xlabel('gl')
axes[0].set_ylabel('gb')
axes[1].set_ylabel('gb')
axes[2].set_ylabel('gb')
axes[3].set_ylabel('gb')
plt.show()


#############################################################
# Finding all Plane points overlapping in GAIS
#############################################################
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

cat = fits.open('starcat_allscans_10-12-18_cuts.fits')[1].data
cat = cat[np.where(cat['expsum'] > 5)]
catgal = SkyCoord(cat['gl']*u.deg, cat['gb']*u.deg, frame='galactic')
cid, gid, angsep, ang3d = gal.search_around_sky(catgal, 1*u.degree)
q = np.unique(cid)
c2 = cat[q]


#############################################################
# SQL query for matching Gaia data from main gaia query page
#############################################################
SELECT crossmatch_positional('user_smohamme','plane','gaiadr2','gaia_source',3.0,'galexmatch') FROM dual

SELECT a."plane_oid", a."col0", a."col1", a."col2", gaia."solution_id", gaia."source_id", gaia."ra", gaia."ra_error", gaia."dec", gaia."dec_error", gaia."parallax", gaia."parallax_error",  gaia."visibility_periods_used", gaia."duplicated_source",gaia."phot_g_mean_flux", gaia."phot_g_mean_flux_error", gaia."phot_g_mean_mag", gaia."phot_bp_mean_flux", gaia."phot_bp_mean_flux_error", gaia."phot_bp_mean_mag", gaia."phot_rp_mean_flux", gaia."phot_rp_mean_flux_error", gaia."phot_rp_mean_mag", gaia."radial_velocity", gaia."radial_velocity_error", gaia."rv_template_teff", gaia."rv_template_logg", gaia."rv_template_fe_h", gaia."phot_variable_flag", gaia."teff_val", gaia."flame_flags", distance(
  POINT('ICRS', a.col1, a.col2),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr2.gaia_source AS gaia, user_smohamme.plane AS a
WHERE 1=CONTAINS(
  POINT('ICRS', a.col1, a.col2),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.0008333333333333334)
)

# Now try for PS1 DR1
SELECT crossmatch_positional('MyDB','planecoords','PanSTARRS_DR1','StackObjectThin',3.0,'galexmatch') FROM dual

SELECT a."ra", a."dec", ps."objID", ps."gra", ps."gdec", ps."gApMag", ps."gApMagErr", ps."rra", ps."rdec", ps."rApMag", ps."rApMagErr",ps."ira", ps."idec", ps."iApMag", ps."iApMagErr",ps."zra", ps."zdec", ps."zApMag", ps."zApMagErr",ps."yra", ps."ydec", ps."yApMag", ps."yApMagErr", distance(
  POINT('ICRS', a.col1, a.col2),
  POINT('ICRS', ps.gra, ps.gdec)) AS dist
FROM PanSTARRS_DR1.StackObjectThin AS ps, MyDB.planecoords AS a
WHERE 1=CONTAINS(
  POINT('ICRS', a.ra, a.dec),
  CIRCLE('ICRS', ps.ra, ps.dec, 0.0008333333333333334)
)

# casjobs says do this except they don't have spgetneighbors...wtf
CREATE TABLE #UPLOAD(
   up_ra FLOAT,
   up_dec FLOAT,
   up_id int
)
INSERT INTO #UPLOAD
SELECT RA AS UP_RA,DEC AS UP_DEC,search_id AS UP_ID
FROM MyDB.planecoords 
CREATE TABLE #tmp (
              up_id int,
               objid bigint
)
INSERT INTO #tmp
EXEC spGetNeighbors 0.05
INSERT INTO MyDB.pan
select a.*,t.objid as matched_id from #tmp t, MyDB.planecoords a  where t.up_id = a.search_id 

#############################################################
# Match dust to plane data. Do twice because query is killed
#############################################################
# 1. Make coord file from plane data from ra and dec
# 2. Upload to gaia archive and set ra and dec in the table
# 3. Cross match using the match button
# 4. Run query from above SQL query section
# 5. Do below twice

# Remember to rename the 'dist' column too 'angsep'
from dustmaps.bayestar import BayestarQuery
cg = fits.open('plane_gaiadr2_06_12_19.fits')[1].data
cg = Table(cg)
cg.rename_column('dist', 'angsep')

#cg.remove_columns(('col0', 'col1', 'col2'))
cg.rename_column('ALPHA_J2000', 'ra_plane')
cg.rename_column('DELTA_J2000', 'dec_plane')
#cg.rename_column('dec', 'dec_gaia')
#cg.rename_column('ra', 'ra_gaia')

negpar = np.where((cg['parallax'] > 0.) & (cg['phot_g_mean_mag'] > 0.) & (cg['phot_bp_mean_mag'] > 0.) & (cg['phot_rp_mean_mag'] > 0.) & (cg['expsum'] > 3.) & (cg['ctsum'] > 1) & (cg['bkgdsum'] > 0) & (cg['angsep']*3600. < 2))
cg = cg[negpar]


cg = cg[1237628:]

pc = 1000./cg['parallax']
distmod = 5. * np.log10(pc) - 5
cggal = SkyCoord(cg['ra_gaia']*u.deg, cg['dec_gaia']*u.deg, distance=pc*u.pc, frame='icrs')

print 'running dust query'
baye = BayestarQuery()
ebv = baye(cggal, mode='median')

cg['dist'] = pc
cg['distmod'] = distmod
cg['ebv'] = ebv
ascii.write(cg, 'plane_gaiadr2_dust_pt2.txt', format='basic')
addhash('plane_gaiadr2_dust_pt2.txt')


# Pickles pysynphot files
# Use site to find out which star to observe. Ex: pickles_uk_7 is a B8V star
# http://www.stsci.edu/hst/observatory/crds/pickles_atlas.html

os.system('export PYSYN_CDBS=~/cdbs/')

filename = os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_27.fits')
sp = S.FileSpectrum(filename)

obs = S.Observation(sp, S.ObsBandpass('galex,nuv'))

# Get ABmag
obs.effstim('abmag')

nuv = []
g = []
r = []
for star in range(1, 132):
    sp = S.FileSpectrum(os.path.join(os.environ['PYSYN_CDBS'], 'grid', 'pickles', 'dat_uvk', 'pickles_uk_'+str(star)+'.fits'))
    nuvobs = S.Observation(sp, S.ObsBandpass('galex,nuv'))
    gobs = S.Observation(sp, S.ObsBandpass('sdss,g'))
    robs = S.Observation(sp, S.ObsBandpass('sdss,r'))

    nuv.append(nuvobs.effstim('abmag'))
    g.append(gobs.effstim('abmag'))
    r.append(robs.effstim('abmag'))

nuv = np.array(nuv)
g = np.array(g)
r = np.array(r)



#############################################################
# Get plane table with no outside catalog matches 
#############################################################
cat = fits.open('starcat_allscans_10-12-18_cuts.fits')[1].data
cg = fits.open('plane_gaiadr2_dust_11_14_18.fits')[1].data

catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
cggal = SkyCoord(cg['ra_plane']*u.deg, cg['dec_plane']*u.deg, frame='icrs')
catind, cgind, angsep, ang3d = search_around_sky(catgal, cggal, 1*u.arcsec)
cat = Table(cat)
cat.remove_rows(catind)

catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
ps = fits.open('plane_ps1_g10-20_10_12_18.fits')[1].data
psgal = SkyCoord(ps['ALPHA_J2000']*u.deg, ps['DELTA_J2000']*u.deg, frame='icrs')
catind, psind, angsep, ang3d = search_around_sky(catgal, psgal, 1*u.arcsec)
cat.remove_rows(catind)

gg = fits.open('plane_gais_11_27_18.fits')[1].data
catgal = SkyCoord(cat['ALPHA_J2000']*u.deg, cat['DELTA_J2000']*u.deg, frame='icrs')
gggal = SkyCoord(gg['ALPHA_J2000']*u.deg, gg['DELTA_J2000']*u.deg, frame='icrs')
catind, ggind, angsep, ang3d = search_around_sky(catgal, gggal, 1*u.arcsec)
cat.remove_rows(catind)
ascii.write(cat, 'plane_innosurveys_11_27_18.txt', format='basic')



# 1. Open exp file, get wcs and convert from pix to wcs
# 2. use scan number to select only single scans in that region
# 3. Convert a list of those scans to pixels
# 4. Get table of expvalues 


from astropy.wcs import wcs
scans = ['0-10', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']

sing = np.array([0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 5.9, 6.8, 8.6, 9.5, 10.4, 11.3, 12.2, 14.0, 14.9, 15.8, 16.7, 17.6, 18.5, 19.4, 20.3, 21.2, 22.1, 23.0, 23.9, 24.8, 25.7, 28.4, 29.3, 30.2, 31.1, 32.0, 32.9, 33.8, 34.7, 35.6, 39.2, 42.8, 43.7, 44.6, 45.5, 46.4, 47.3, 48.2, 49.1, 50.0, 67.1, 68.9, 71.6, 74.3, 75.2, 76.1, 77.0, 77.9, 78.8, 79.7, 80.6, 81.5, 82.4, 83.3, 87.8, 88.7, 89.6, 90.5, 91.4, 92.3, 93.2, 94.1, 95.0, 95.9, 96.8, 97.7, 98.6, 99.5, 100.4, 101.3, 102.2, 103.1, 104.0, 104.9, 105.8, 106.7, 107.6, 110.3, 111.2, 112.1, 113.0, 113.9, 114.8, 119.3, 121.1, 122.9, 124.7, 125.6, 126.5, 127.4, 128.3, 129.2, 130.1, 131.0, 131.9, 132.8, 133.7, 134.6, 135.5, 136.4, 137.3, 138.2, 139.1, 140.0, 140.9, 141.8, 143.6, 144.5, 145.4, 148.1, 149.0, 149.9, 150.8, 151.7, 152.6, 153.5, 155.3, 156.2, 157.1, 158.0, 160.7, 161.6, 163.4, 167.0, 167.9, 172.4, 173.3, 174.2, 175.1, 176.0, 176.9, 177.8, 178.7, 179.6, 180.5, 183.2, 185.0, 190.4, 191.3, 197.6, 198.5, 200.3, 201.2, 203.0, 203.9, 205.7, 206.6, 207.5, 208.4, 209.3, 210.2, 211.1, 212.0, 212.9, 213.8, 214.7, 215.6, 216.5, 217.4, 218.3, 219.2, 220.1, 221.0, 221.9, 222.8, 223.7, 224.6, 225.5, 226.4, 228.2, 229.1, 230.0, 230.9, 231.8, 234.5, 235.4, 236.3, 237.2, 238.1, 239.0, 239.9, 240.8, 241.7, 242.6, 243.5, 244.4, 245.3, 246.2, 247.1, 248.0, 248.9, 249.8, 250.7, 251.6, 252.5, 253.4, 254.3, 255.2, 256.1, 257.0, 258.8, 259.7, 260.6, 261.5, 263.3, 264.2, 265.1, 266.0, 266.9, 268.7, 269.6, 270.5, 271.4, 272.3, 273.2, 274.1, 275.0, 275.9, 276.8, 278.6, 279.5, 281.3, 283.1, 284.0, 285.8, 286.7, 288.5, 289.4, 290.3, 291.2, 292.1, 293.0, 293.9, 295.7, 297.5, 298.4, 301.1, 302.0, 302.9, 303.8, 304.7, 305.6, 306.5, 308.3, 309.2, 310.1, 315.5, 316.4, 317.3, 318.2, 319.1, 320.0, 320.9, 321.8, 322.7, 323.6, 324.5, 325.4, 326.3, 327.2, 328.1, 329.0, 329.9, 331.7, 332.6, 333.5, 334.4, 335.3, 338.0, 338.9, 339.8, 341.6, 342.5, 343.4, 345.2, 348.8, 349.7, 350.6, 351.5, 352.4, 353.3, 354.2, 355.1, 356.0, 357.8, 358.7, 359.0])


allexp = Table()

for scan in scans:
    print scan
    scanexp = Table()
    exp = fits.open('../galexscans/count_map_'+scan+'_exp.fits')[0].data

    hdulist = fits.open('../galexscans/count_map_'+scan+'_exp.fits')
    w = wcs.WCS(hdulist[0].header)
    pixels = np.array([[0, np.shape(exp)[1]], [0, np.shape(exp)[0]]]).T
    world = w.wcs_pix2world(pixels, 1)
    
    gllimits = [world[1][0], world[0][0]]
    if scan == '0-10':
        gllimits[0] = 0
    elif scan == '351-1':
        gllimits[1] = 360
    glrange = np.linspace(gllimits[0], gllimits[1], (np.shape(exp)[0]))
    
    region = sing[np.where((sing > gllimits[0]) & (sing < gllimits[1]))]
    pixvals = w.wcs_world2pix(region, np.zeros(len(region)), 1)[0]
    
    for line in range(len(pixvals)):
        a = exp[:, int(pixvals[line])]
        scanexp['scan'+str(region[line])] = a
    
    #ascii.write(scanexp, 'expdata_'+scan+'.txt', format='basic')
    #addhash('expdata_'+scan+'.txt')
    allexp = hstack([allexp, scanexp])
ascii.write(allexp, 'expdata_11_28_18.txt', format='basic')


exp = fits.open('expdata_11_28_18.fits')[1].data

scans = ['scan0.5', 'scan1.4', 'scan2.3', 'scan3.2', 'scan4.1', 'scan5.0', 'scan5.9', 'scan6.8', 'scan8.6', 'scan9.5', 'scan10.4', 'scan11.3', 'scan12.2', 'scan14.0', 'scan14.9', 'scan15.8', 'scan16.7', 'scan17.6', 'scan18.5', 'scan19.4', 'scan20.3', 'scan21.2', 'scan22.1', 'scan23.0', 'scan23.9', 'scan24.8', 'scan25.7', 'scan28.4', 'scan29.3', 'scan30.2', 'scan31.1', 'scan32.0', 'scan32.9', 'scan33.8', 'scan34.7', 'scan35.6', 'scan39.2', 'scan42.8', 'scan43.7', 'scan44.6', 'scan45.5', 'scan46.4', 'scan47.3', 'scan48.2', 'scan49.1', 'scan50.0', 'scan67.1', 'scan68.9', 'scan71.6', 'scan74.3', 'scan75.2', 'scan76.1', 'scan77.0', 'scan77.9', 'scan78.8', 'scan79.7', 'scan80.6', 'scan81.5', 'scan82.4', 'scan83.3', 'scan87.8', 'scan88.7', 'scan89.6', 'scan90.5', 'scan91.4', 'scan92.3', 'scan93.2', 'scan94.1', 'scan95.0', 'scan95.9', 'scan96.8', 'scan97.7', 'scan98.6', 'scan99.5', 'scan100.4', 'scan101.3', 'scan102.2', 'scan103.1', 'scan104.0', 'scan104.9', 'scan105.8', 'scan106.7', 'scan107.6', 'scan110.3', 'scan111.2', 'scan112.1', 'scan113.0', 'scan113.9', 'scan114.8', 'scan119.3', 'scan121.1', 'scan122.9', 'scan124.7', 'scan125.6', 'scan126.5', 'scan127.4', 'scan128.3', 'scan129.2', 'scan130.1', 'scan131.0', 'scan131.9', 'scan132.8', 'scan133.7', 'scan134.6', 'scan135.5', 'scan136.4', 'scan137.3', 'scan138.2', 'scan139.1', 'scan140.0', 'scan140.9', 'scan141.8', 'scan143.6', 'scan144.5', 'scan145.4', 'scan148.1', 'scan149.0', 'scan149.9', 'scan150.8', 'scan151.7', 'scan152.6', 'scan153.5', 'scan155.3', 'scan156.2', 'scan157.1', 'scan158.0', 'scan160.7', 'scan161.6', 'scan163.4', 'scan167.0', 'scan167.9', 'scan172.4', 'scan173.3', 'scan174.2', 'scan175.1', 'scan176.0', 'scan176.9', 'scan177.8', 'scan178.7', 'scan179.6', 'scan180.5', 'scan183.2', 'scan185.0', 'scan190.4', 'scan191.3', 'scan197.6', 'scan198.5', 'scan200.3', 'scan201.2', 'scan203.0', 'scan203.9', 'scan205.7', 'scan206.6', 'scan207.5', 'scan208.4', 'scan209.3', 'scan210.2', 'scan211.1', 'scan212.0', 'scan212.9', 'scan213.8', 'scan214.7', 'scan215.6', 'scan216.5', 'scan217.4', 'scan218.3', 'scan219.2', 'scan220.1', 'scan221.0', 'scan221.9', 'scan222.8', 'scan223.7', 'scan224.6', 'scan225.5', 'scan226.4', 'scan228.2', 'scan229.1', 'scan230.0', 'scan230.9', 'scan231.8', 'scan234.5', 'scan235.4', 'scan236.3', 'scan237.2', 'scan238.1', 'scan239.0', 'scan239.9', 'scan240.8', 'scan241.7', 'scan242.6', 'scan243.5', 'scan244.4', 'scan245.3', 'scan246.2', 'scan247.1', 'scan248.0', 'scan248.9', 'scan249.8', 'scan250.7', 'scan251.6', 'scan252.5', 'scan253.4', 'scan254.3', 'scan255.2', 'scan256.1', 'scan257.0', 'scan258.8', 'scan259.7', 'scan260.6', 'scan261.5', 'scan263.3', 'scan264.2', 'scan265.1', 'scan266.0', 'scan266.9', 'scan268.7', 'scan269.6', 'scan270.5', 'scan271.4', 'scan272.3', 'scan273.2', 'scan274.1', 'scan275.0', 'scan275.9', 'scan276.8', 'scan278.6', 'scan279.5', 'scan281.3', 'scan283.1', 'scan284.0', 'scan285.8', 'scan286.7', 'scan288.5', 'scan289.4', 'scan290.3', 'scan291.2', 'scan292.1', 'scan293.0', 'scan293.9', 'scan295.7', 'scan297.5', 'scan298.4', 'scan301.1', 'scan302.0', 'scan302.9', 'scan303.8', 'scan304.7', 'scan305.6', 'scan306.5', 'scan308.3', 'scan309.2', 'scan310.1', 'scan315.5', 'scan316.4', 'scan317.3', 'scan318.2', 'scan319.1', 'scan320.0', 'scan320.9', 'scan321.8', 'scan322.7', 'scan323.6', 'scan324.5', 'scan325.4', 'scan326.3', 'scan327.2', 'scan328.1', 'scan329.0', 'scan329.9', 'scan331.7', 'scan332.6', 'scan333.5', 'scan334.4', 'scan335.3', 'scan338.0', 'scan338.9', 'scan339.8', 'scan341.6', 'scan342.5', 'scan343.4', 'scan345.2', 'scan348.8', 'scan349.7', 'scan350.6', 'scan351.5', 'scan352.4', 'scan353.3', 'scan354.2', 'scan355.1', 'scan356.0', 'scan357.8', 'scan358.7', 'scan359.0']

expvals = 0
for i in range(len(scans)):
    expvals += exp[scans[i]]
    #if (i % 20 == 0) & (i != 0):
    expvals = expvals#/20.
    plt.plot(np.linspace(-10, 10, len(expvals)), expvals, label=str(i), linewidth=0.1, alpha=0.3)
    expvals = 0
plt.xlabel('Galactic Latitude')
plt.ylabel('Exposure time [s]')
plt.ylim(0, 200)
plt.show()



plt.ylim(14, -3)
plt.xlim(-2, 11.6)
plt.xlabel('(NUV - G)$_0$')
plt.ylabel('M$_G$')




psg.remove_columns(('NUMBER_2','X_IMAGE_2','Y_IMAGE_2','ra_plane','dec_plane','FLUX_AUTO_2','FLUXERR_AUTO_2','FLUX_APER_2','A_IMAGE_2','B_IMAGE_2','THETA_IMAGE_2','FWHM_IMAGE_2','nuv_2','gl_2','gb_2','expsum_2','ctsum_2','bkgdsum_2'))



cgp = fits.open('plane_gaia_ps2_03-22-19.fits')[1].data
nuvg = (cgp['nuv']-cgp['ebv']*7.24)-(cgp['phot_g_mean_mag']-cgp['ebv']*2.85)
mg = cgp['phot_g_mean_mag']-cgp['distmod']-cgp['ebv']*2.85
obcut = np.where((mg < -0.5) & (mg > -2) & (nuvg > 0) & (nuvg < 1))
ob = cgp[obcut]
cuts,= np.where((ob['visibility_periods_used'] > 8) & (ob['expsum'] > 10) & (ob['parallax_error']/ob['parallax'] < 0.1) & (ob['ng'] > 8) & (ob['nr'] > 8) & (ob['ni'] > 8) & (ob['nz'] > 8) & (ob['ny'] > 8))
ob2 = ob[cuts]



q = np.where((cg['gl'] > mol['llower'][n]) & (cg['gl'] < mol['lupper'][n]) & (cg['gb'] > mol['blower'][n]) & (cg['gb'] > mol['bupper'][n]) & (cg['dist'] > mol['dist'][n]-100) & (cg['dist'] < mol['dist'][n]+100))


plt.hist(c1['nuv'], range=[12,20], bins=20, label=mol['name'][1])
plt.hist(c0['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][0])
plt.hist(c5['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][5])
plt.hist(c10['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][10])
plt.hist(c11['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][11])
plt.hist(c4['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][4])
plt.hist(c7['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][7])
plt.hist(c9['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][9])
plt.hist(c12['nuv'], range=[12,20], bins=20, histtype='step', stacked=True, label=mol['name'][12])
plt.legend(loc=2)
plt.xlabel('NUV')
plt.show()




plt.hist((c1['nuv']-c1['ebv']*7.24)-(c1['phot_g_mean_mag']-c1['ebv']*2.85), range=[0, 10], bins=20, label=mol['name'][1])
plt.hist((c0['nuv']-c0['ebv']*7.24)-(c0['phot_g_mean_mag']-c0['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][0])
plt.hist((c4['nuv']-c4['ebv']*7.24)-(c4['phot_g_mean_mag']-c4['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][4])
plt.hist((c5['nuv']-c5['ebv']*7.24)-(c5['phot_g_mean_mag']-c5['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][5])
plt.hist((c7['nuv']-c7['ebv']*7.24)-(c7['phot_g_mean_mag']-c7['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][7])
plt.hist((c9['nuv']-c9['ebv']*7.24)-(c9['phot_g_mean_mag']-c9['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][9])
plt.hist((c10['nuv']-c10['ebv']*7.24)-(c10['phot_g_mean_mag']-c10['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][10])
plt.hist((c11['nuv']-c11['ebv']*7.24)-(c11['phot_g_mean_mag']-c11['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][11])
plt.hist((c12['nuv']-c12['ebv']*7.24)-(c12['phot_g_mean_mag']-c12['ebv']*2.85), range=[0, 10], bins=20, histtype='step', stacked=True, label=mol['name'][12])
plt.legend(loc=2)
plt.xlabel('(NUV - G)$_0$')
plt.show()



cg = fits.open('plane_gaiadr2_dust_06_12_19.fits')[1].data
cg = cg[np.where((cg['expsum'] > 3) & (cg['ebv'] > 0))]
    
plt.hist2d(cg['nuv']-cg['phot_g_mean_mag'], cg['phot_g_mean_mag']-cg['distmod'], bins=(1000,1000), norm=matplotlib.colors.LogNorm())
plt.xlim(-2, 11.6)
plt.ylim(14, -3)
plt.show()

nuvg = cg['nuv'] - cg['phot_g_mean_mag']
bprp = cg['phot_bp_mean_mag'] - cg['phot_rp_mean_mag']
mg = cg['phot_g_mean_mag'] - cg['distmod']
rc1cut = np.where((bprp > 1.11) & (mg < 1.11))
rc2cut = np.where((nuvg > 6.8) & (mg < 1.4))





