#######################################################################################
# CMD
#######################################################################################
from mpl_toolkits.axes_grid.inset_locator import inset_axes

ga = Table.read('gais_gaiadr2_apogee_dust_7-29.txt', format='ascii')
pa = Table.read('plane_gaiadr2_dust_apogee_7-28.txt', format='ascii')

lg = Table.read('gais_gaia_lamost_huang20.txt', format='ascii')
lp = Table.read('rc_plane_gaiadr2_lamost_08-28-20.txt', format='ascii')



fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
fig.subplots_adjust(wspace=0)

ax1.scatter((ga['mag_nuv']-ga['ebvsfd']*7.24)-(ga['phot_g_mean_mag']-ga['ebvsfd']*2.85), (ga['phot_g_mean_mag']-ga['ebvsfd']*2.85)-ga['distmod'], s=1, alpha=0.5)
ax2.scatter((pa['nuv']-pa['ebv']*7.24)-(pa['phot_g_mean_mag']-pa['ebv']*2.85), (pa['phot_g_mean_mag']-pa['ebv']*2.85)-pa['distmod'], s=1, alpha=0.5)
ax1.set_ylim(8, -4)
ax1.set_xlim(0, 11)
plt.ion()
ax1.set_xlabel('(NUV-G)$_0$')
ax2.set_xlabel('(NUV-G)$_0$')
ax1.set_ylabel('M$_{G_0}$')
ax1.annotate('GAIS+APOGEE', xy=(6.6, -0.5))
ax2.annotate('UVGAPS+APOGEE', xy=(6.4, -0.5))
ax1.plot((6.5, 10), (-0.3,-0.3), color='red')
ax1.plot((6.5, 10), (1,1), color='red')
ax1.plot((10, 10), (-0.3, 1), color='red')
ax1.plot((6.5, 6.5), (-0.3, 1), color='red')
ax2.plot((6.5, 10), (-0.3,-0.3), color='red')
ax2.plot((6.5, 10), (1,1), color='red')
ax2.plot((10, 10), (-0.3, 1), color='red')
ax2.plot((6.5, 6.5), (-0.3, 1), color='red')



#axin1 = inset_axes(ax1, height=2, width=2, loc=3)

axin1 = inset_axes(ax1, width="100%", height="100%",
                   bbox_to_anchor=(6/11., 9/10., 0.37, -.12),
                   bbox_transform=ax1.transAxes, loc=2, borderpad=0)

axin1.set_xlim(6, 10)
axin1.set_ylim(1.5, -1.5)
axin1.scatter((lg['mag_nuv']-lg['ebv']*7.24)-(lg['phot_g_mean_mag']-lg['ebv']*2.85), lg['phot_g_mean_mag']-lg['ebv']*2.85-lg['distmod'], s=1, alpha=0.1)
ax1.annotate('GAIS+LAMOST', xy=(6.6, -3))
axin1.invert_yaxis()


axin2 = inset_axes(ax2, width="100%", height="100%",
                   bbox_to_anchor=(6/11., 9/10., 0.37, -.12),
                   bbox_transform=ax2.transAxes, loc=2, borderpad=0)

axin2.set_xlim(6, 10)
axin2.set_ylim(1.5, -1.5)
axin2.scatter((lp['nuv']-lp['ebv3d']*7.24)-(lp['phot_g_mean_mag']-lp['ebv3d']*2.85), lp['phot_g_mean_mag']-lp['ebv3d']*2.85-lp['distmod'], s=1, alpha=0.5)
ax2.annotate('UVGAPS+LAMOST', xy=(6.4, -3))
axin2.invert_yaxis()


#ax2.arrow(4, 5.5, (2.972-0.789)*-1, -0.789, head_length=0.2, head_width=0.2, color='red')
ax2.arrow(4, 7, (7.24-2.85)*-1*.5, -2.85*.5, head_length=0.2, head_width=0.2, color='red')

plt.show()


#######################################################################################
# CMD no dust
#######################################################################################
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
fig.subplots_adjust(wspace=0)
ax1.scatter((ga['mag_nuv'])-(ga['phot_g_mean_mag']), (ga['phot_g_mean_mag'])-ga['distmod'], s=1, alpha=0.5)
ax2.scatter((pa['nuv']-pa['phot_g_mean_mag']), (pa['phot_g_mean_mag']-pa['distmod']), s=1, alpha=0.5)
ax1.set_ylim(10, -2)
ax1.set_xlim(0, 11)
plt.ion()
ax1.set_xlabel('(NUV-G)$_0$')
ax2.set_xlabel('(NUV-G)$_0$')
ax1.set_ylabel('M$_{G_0}$')
ax1.annotate('GAIS+Apogee', xy=(0, 9.9))
ax2.annotate('UVGAPS+Apogee', xy=(0, 9.9))
plt.show()


#######################################################################################
# histograms
#######################################################################################
#m = (0.265-0.065)/(-0.88-0.02)
#b = 0.0694

x1, y1 = 0.02, 0.065
x2, y2 = -0.86, 0.23
m = (y2 - y1) / (x2 - x1)
b = y1 - m*x1

thickp,= np.where((rcp['alphafe'] > 0.08) & (rcp['alphafe'] > (m*rcp['FE_H'] + b)))
thinp,= np.where((rcp['alphafe'] < 0.08) | (rcp['AFE'] < (m*rcp['FE_H'] + b)))
thickg,= np.where((rcg['alphafe'] > 0.08) & (rcg['alphafe'] > (m*rcg['FE_H'] + b)))
thing,= np.where((rcg['alphafe'] < 0.08) | (rcg['alphafe'] < (m*rcg['FE_H'] + b)))
heap beer lyricsthinp = rcp[thinp]
thickp = rcp[thickp]
thing = rcg[thing]
thickg = rcg[thickg]

fig, axes = plt.subplots(2, 4)
axes[0, 0].hist(rcg['glat'], range=[-90,90])
axes[0, 0].hist(thing['glat'], range=[-90,90], histtype='step',color='red', fill=False, stacked=True)
axes[0, 0].hist(thickg['glat'], range=[-90,90], histtype='step', color='black', fill=False, stacked=True)
axes[1, 0].hist(rcp['gb'], range=[-10,10])
axes[1, 0].hist(thinp['gb'], range=[-10,10], histtype='step',color='red', fill=False, stacked=True)
axes[1, 0].hist(thickp['gb'], range=[-10,10], histtype='step', color='black', fill=False, stacked=True)
axes[0, 1].hist(np.log10(rcg['dist']), range=[1.8,3.5], label='All RC')
axes[0, 1].hist(np.log10(thing['dist']), range=[1.8,3.5], histtype='step', color='red', stacked=True, fill=False, label='Thin Disk')
axes[0, 1].hist(np.log10(thickg['dist']), range=[1.8,3.5], histtype='step', color='black', stacked=True, fill=False, label='Thick Disk')

axes[1, 1].hist(np.log10(rcp['dist']), range=[1.8,3.5])
axes[1, 1].hist(np.log10(thinp['dist']), range=[1.8,3.5], histtype='step', color='red', stacked=True, fill=False, label='Thin disk')
axes[1, 1].hist(np.log10(thickp['dist']), range=[1.8,3.5], histtype='step', color='black', stacked=True, fill=False, label='Thick disk')
axes[0, 2].hist(rcg['mag_nuv'], range=[12.5, 22.5])
axes[0, 2].hist(thing['mag_nuv'],range=[12.5, 22.5], histtype='step', color='red', stacked=True, fill=False)
axes[0, 2].hist(thickg['mag_nuv'],range=[12.5, 22.5], histtype='step', color='black', stacked=True, fill=False)

axes[1, 2].hist(rcp['nuv'], range=[12.5, 22.5])
axes[1, 2].hist(thinp['nuv'],range=[12.5, 22.5], histtype='step', color='red', stacked=True, fill=False)
axes[1, 2].hist(thickp['nuv'],range=[12.5, 22.5], histtype='step', color='black', stacked=True, fill=False)
axes[0, 3].hist(np.log10(rcg['ebvsfd']), range=[-3,0])
axes[0, 3].hist(np.log10(thing['ebvsfd']), range=[-3,0], histtype='step', color='red', stacked=True, fill=False)
axes[0, 3].hist(np.log10(thickg['ebvsfd']), range=[-3,0], histtype='step', color='black', stacked=True, fill=False)

axes[1, 3].hist(np.log10(rcp['ebv3d']), range=[-3,0])
axes[1, 3].hist(np.log10(thinp['ebv3d']), range=[-3,0], histtype='step', color='red', stacked=True, fill=False)
axes[1, 3].hist(np.log10(thickp['ebv3d']), range=[-3,0], histtype='step', color='black', stacked=True, fill=False)
axes[1, 0].set_xlabel('Galactic Latitude')
axes[1, 1].set_xlabel('log Distance [pc]')
axes[1, 2].set_xlabel('NUV')
axes[1, 3].set_xlabel('log E(B-V)')
axes[0, 1].legend()
axes[0, 0].annotate('GAIS', xy=(-95, 4500))
axes[1, 0].annotate('UVGAPS', xy=(-10, 140))
#fig.subplots_adjust(hspace=0)
plt.rcParams.update({'font.size': 14})
plt.show()




######################################################################
# Color vs TEFF vs Fe/H
######################################################################
#m = (0.265-0.065)/(-0.88-0.02)
#b = 0.0694


x1, y1 = 0.02, 0.065
x2, y2 = -0.86, 0.23
m = (y2 - y1) / (x2 - x1)
b = y1 - m*x1


thickp,= np.where((rcp['alphafe'] > 0.08) & (rcp['alphafe'] > (m*rcp['FE_H'] + b)))
thinp,= np.where((rcp['alphafe'] < 0.08) | (rcp['alphafe'] < (m*rcp['FE_H'] + b)))
thickg,= np.where((rcg['alphafe'] > 0.08) & (rcg['alphafe'] > (m*rcg['FE_H'] + b)))
thing,= np.where((rcg['alphafe'] < 0.08) | (rcg['alphafe'] < (m*rcg['FE_H'] + b)))
thinp = rcp[thinp]
thickp = rcp[thickp]
thing = rcg[thing]
thickg = rcg[thickg]


#y1 = (rc['B_apass']-rc['ebv']*3.626)-(rc['V_apass']-rc['ebv']*2.742)
y1 = rc['phot_bp_mean_mag']-rc['phot_rp_mean_mag']
y2 = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)
x = rc['TEFF']
m1, b1, rval1, pval1, stderr1 = stats.linregress(x, y1)
line1 = m1*x + b1
err1 = np.sqrt(np.sum((line1-y1)**2/len(y1)))
m2, b2, rval2, pval2, stderr2 = stats.linregress(x, y2)
line2 = m2*x + b2
err2 = np.sqrt(np.sum((line2-y2)**2/len(y2)))

bprperr = np.sqrt(rc['BPerr']**2+rc['RPerr']**2)
nuvgerr = np.sqrt(rc['nuv_magerr']**2+rc['Gerr']**2)

fig, (ax1, ax2) = plt.subplots(2, 1)
cmap = ax1.scatter(x[thin], y1[thin], s=40, c=rc['FE_H'][thin], cmap=cm.viridis_r, vmin=-.5, vmax=.35, label='Thin disk', edgecolor='black')
ax1.scatter(x[thick], y1[thick], s=40, c=rc['FE_H'][thick], cmap=cm.viridis_r, vmin=-.5, vmax=.35, label='Thick disk', edgecolor='black', linewidth=1)
ax1.errorbar(x, y1, xerr=rc['TEFF_ERR'], yerr=bprperr, ecolor='gray', fmt='none', marker='none', mew=0, elinewidth=1.3, **{"zorder":0})
ax1.plot(x, line1, linewidth=3, c='red', zorder=10)

ax2.scatter(x[thin], y2[thin], s=40, c=rc['FE_H'][thin], cmap=cm.viridis_r, vmin=-.5, vmax=.35, edgecolor='black')
ax2.scatter(x[thick], y2[thick], s=40, c=rc['FE_H'][thick], cmap=cm.viridis_r, vmin=-.5, vmax=.35, edgecolor='black', linewidth=1)
ax2.errorbar(x, y2, xerr=rc['TEFF_ERR'], yerr=nuvgerr, ecolor='gray', fmt='none', marker='none', mew=0, elinewidth=1.3, **{"zorder":0})
ax2.plot(x, line2, linewidth=3, c='red', zorder=10)

ax1.set_ylabel('(G$_{BP}$ - G$_{RP}$)$_0$')
ax2.set_xlabel('T$_{eff}$ [K]')
ax2.set_ylabel('(NUV - G)$_0$')
fig.subplots_adjust(right=.84)
fig.subplots_adjust(hspace=0)
ax1.set_xticklabels([])

ax1.annotate('(G$_{BP}$ - G$_{RP}$)$_0$ = -0.000362 * T$_{eff}$ + 3.0', xy=(5078, 0.55), color='black', size=13)
ax1.annotate('$\sigma$ = '+str(round(err1, 3)), xy=(5253, 0.62), color='black', size=13)

ax2.annotate('(NUV - G)$_0$ = -0.00474 * T$_{eff}$ + 30.958', xy=(5078, 3.2), color='black', size=13)
ax2.annotate('$\sigma$ = '+str(round(err2, 3)), xy=(5253, 3.5), color='black', size=13)

cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
plt.show()


#######################################################################################
# Teff vs bprp and nuvg, plane vs gais
#######################################################################################
xp = np.linspace(3800, 5300, 50)

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thickp,= np.where((rcp['alphafe'] > 0.08) & (rcp['alphafe'] > (m*rcp['FE_H'] + b)))
thinp,= np.where((rcp['alphafe'] < 0.08) | (rcp['alphafe'] < (m*rcp['FE_H'] + b)))
thickg,= np.where((rcg['alphafe'] > 0.08) & (rcg['alphafe'] > (m*rcg['FE_H'] + b)))
thing,= np.where((rcg['alphafe'] < 0.08) | (rcg['alphafe'] < (m*rcg['FE_H'] + b)))
thinp = rcp[thinp]
thickp = rcp[thickp]
thing = rcg[thing]
thickg = rcg[thickg]



fig, axes = plt.subplots(2, 2)
axes[0, 0].scatter(rcg['TEFF'], rcg['phot_bp_mean_mag'] - rcg['phot_rp_mean_mag'], c=rcg['FE_H'], vmin=-0.5, vmax=0.4, s=4)
axes[0, 1].scatter(rcp['TEFF'], rcp['phot_bp_mean_mag'] - rcp['phot_rp_mean_mag'], c=rcp['FE_H'], vmin=-0.5, vmax=0.4, s=4)
axes[0, 0].plot(xp, -0.000362*xp + 3, c='red')
axes[0, 1].plot(xp, -0.000362*xp + 3, c='red')
axes[1, 0].scatter(rcg['TEFF'], (rcg['mag_nuv']-rcg['ebvsfd']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebvsfd']*2.85), c=rcg['FE_H'], vmin=-0.5, vmax=0.4, s=4)
cmap = axes[1, 1].scatter(rcp['TEFF'], (rcp['nuv']-rcp['ebv']*7.24)-(rcp['phot_g_mean_mag']-rcp['ebv']*2.85), c=rcp['FE_H'], vmin=-0.5, vmax=0.4, s=4)


cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')

axes[1, 0].plot(xp, -0.00474*xp + 30.958, c='red')
axes[1, 1].plot(xp, -0.00474*xp + 30.958, c='red')

axes[1, 0].set_xlabel('T$_{eff}$')
axes[1, 1].set_xlabel('T$_{eff}$')
axes[0, 0].set_ylabel('(BP-RP)$_0$')
axes[1, 0].set_ylabel('(NUV-G)$_0$')


axes[0, 0].set_xlim(4200, 5400)
axes[0, 1].set_xlim(4200, 5400)
axes[1, 0].set_xlim(4200, 5400)
axes[1, 1].set_xlim(4200, 5400)

axes[0, 0].set_ylim(0.9, 1.8)
axes[0, 1].set_ylim(0.9, 1.8)
axes[1, 0].set_ylim(5, 10.5)
axes[1, 1].set_ylim(5, 10.5)
fig.subplots_adjust(hspace=0)

axes[0, 0].annotate('GAIS', xy=(5000, 2))
axes[0, 1].annotate('UVGAPS', xy=(5000, 2))

plt.show()


#######################################################################################
# Teff vs nuvg, plane vs gais, thin and thick
#######################################################################################
#m = (0.265-0.065)/(-0.88-0.02)
#b = 0.0694

rcp = Table.read('rc_plane_gaiadr2_apogee_7-28.txt', format='ascii')
rcg = Table.read('rc_gais_gaiadr2_apogee_dust_7-29.txt', format='ascii')

x1, y1 = 0.02, 0.065
x2, y2 = -0.86, 0.23
m = (y2 - y1) / (x2 - x1)
b = y1 - m*x1

thickp,= np.where((rcp['alphafe'] > 0.08) & (rcp['alphafe'] > (m*rcp['FE_H'] + b)))
thinp,= np.where((rcp['alphafe'] < 0.08) | (rcp['alphafe'] < (m*rcp['FE_H'] + b)))
thickg,= np.where((rcg['alphafe'] > 0.08) & (rcg['alphafe'] > (m*rcg['FE_H'] + b)))
thing,= np.where((rcg['alphafe'] < 0.08) | (rcg['alphafe'] < (m*rcg['FE_H'] + b)))
thinp = rcp[thinp]
thickp = rcp[thickp]
thing = rcg[thing]
thickg = rcg[thickg]
xline = np.linspace(3800, 5300, 50)



# Get fit lines for plot
xg = rcg['TEFF']
xp = rcp['TEFF']
xgthick = thickg['TEFF']
xpthick = thickp['TEFF']
xgthin = thing['TEFF']
xpthin = thinp['TEFF']


yg = (rcg['mag_nuv']-rcg['ebv3d']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebv3d']*2.85)
yp = (rcp['nuv']-rcp['ebv']*7.24)-(rcp['phot_g_mean_mag']-rcp['ebv']*2.85)
ygthick = (thickg['mag_nuv']-thickg['ebv3d']*7.24)-(thickg['phot_g_mean_mag']-thickg['ebv3d']*2.85)
ypthick = (thickp['nuv']-thickp['ebv']*7.24)-(thickp['phot_g_mean_mag']-thickp['ebv']*2.85)
ygthin = (thing['mag_nuv']-thing['ebv3d']*7.24)-(thing['phot_g_mean_mag']-thing['ebv3d']*2.85)
ypthin = (thinp['nuv']-thinp['ebv']*7.24)-(thinp['phot_g_mean_mag']-thinp['ebv']*2.85)

zg, vg = np.polyfit(xg, yg, 1, cov=True)
zp, vp = np.polyfit(xp, yp, 1, cov=True)
pg = np.poly1d(zg)
pp = np.poly1d(zp)

zgthick, vgthick = np.polyfit(xgthick, ygthick, 1, cov=True)
zpthick, vpthick = np.polyfit(xpthick, ypthick, 1, cov=True)
pgthick = np.poly1d(zgthick)
ppthick = np.poly1d(zpthick)

zgthin, vgthin = np.polyfit(xgthin, ygthin, 1, cov=True)
zpthin, vpthin = np.polyfit(xpthin, ypthin, 1, cov=True)
pgthin = np.poly1d(zgthin)
ppthin = np.poly1d(zpthin)

gerr = np.sqrt(np.sum((pg(xg)-yg)**2)/len(xg))
perr = np.sqrt(np.sum((pp(xp)-yp)**2)/len(xp))

gthickerr = np.sqrt(np.sum((pgthick(xgthick)-ygthick)**2)/len(xgthick))
pthickerr = np.sqrt(np.sum((ppthick(xpthick)-ypthick)**2)/len(xpthick))

gthinerr = np.sqrt(np.sum((pgthin(xgthin)-ygthin)**2)/len(xgthin))
pthinerr = np.sqrt(np.sum((ppthin(xpthin)-ypthin)**2)/len(xpthin))




fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
axes[0, 0].scatter(rcg['TEFF'], (rcg['mag_nuv']-rcg['ebv3d']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebv3d']*2.85), c=rcg['FE_H'], vmin=-0.5, vmax=0.4, s=4)
cmap = axes[0, 1].scatter(rcp['TEFF'], (rcp['nuv']-rcp['ebv']*7.24)-(rcp['phot_g_mean_mag']-rcp['ebv']*2.85), c=rcp['FE_H'], vmin=-0.5, vmax=0.4, s=20)

axes[1, 0].scatter(thing['TEFF'], (thing['mag_nuv']-thing['ebv3d']*7.24)-(thing['phot_g_mean_mag']-thing['ebv3d']*2.85), c=thing['FE_H'], vmin=-0.5, vmax=0.4, s=4)
axes[1, 1].scatter(thinp['TEFF'], (thinp['nuv']-thinp['ebv']*7.24)-(thinp['phot_g_mean_mag']-thinp['ebv']*2.85), c=thinp['FE_H'], vmin=-0.5, vmax=0.4, s=20)

axes[2, 0].scatter(thickg['TEFF'], (thickg['mag_nuv']-thickg['ebv3d']*7.24)-(thickg['phot_g_mean_mag']-thickg['ebv3d']*2.85), c=thickg['FE_H'], vmin=-0.5, vmax=0.4, s=4)
axes[2, 1].scatter(thickp['TEFF'], (thickp['nuv']-thickp['ebv']*7.24)-(thickp['phot_g_mean_mag']-thickp['ebv']*2.85), c=thickp['FE_H'], vmin=-0.5, vmax=0.4, s=20)

axes[0, 0].plot(xp, -0.00474*xp + 30.958, c='red')
axes[0, 1].plot(xp, -0.00474*xp + 30.958, c='red')
axes[1, 0].plot(xp, -0.00474*xp + 30.958, c='red')
axes[1, 1].plot(xp, -0.00474*xp + 30.958, c='red')
axes[2, 0].plot(xp, -0.00474*xp + 30.958, c='red')
axes[2, 1].plot(xp, -0.00474*xp + 30.958, c='red')
cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
axes[2, 0].set_xlabel('T$_{eff}$')
axes[2, 1].set_xlabel('T$_{eff}$')
axes[0, 0].set_ylabel('(NUV-G)$_0$')
axes[1, 0].set_ylabel('(NUV-G)$_0$')
axes[2, 0].set_ylabel('(NUV-G)$_0$')
axes[0, 0].set_xlim(4200, 5400)
axes[0, 1].set_xlim(4200, 5400)
axes[1, 0].set_xlim(4200, 5400)
axes[1, 1].set_xlim(4200, 5400)
axes[1, 0].set_ylim(6, 10.5)
axes[1, 1].set_ylim(6, 10.5)
fig.subplots_adjust(wspace=0, hspace=0)
axes[0, 0].annotate('GAIS', xy=(4250, 6.5))
axes[0, 1].annotate('UVGAPS', xy=(4250, 6.5))
axes[1, 0].annotate('Thin', xy=(4250, 6.5))
axes[2, 0].annotate('Thick', xy=(4250, 6.5))

axes[0, 0].annotate('(NUV-G)$_0$ ='+str(zg[0])[:7]+'*T$_{eff}$ '+str(zg[1])[:4], xy=(5000, 10))
axes[0, 1].annotate('(NUV-G)$_0$ ='+str(zp[0])[:7]+'*T$_{eff}$ '+str(zp[1])[:4], xy=(5000, 10))
axes[1, 0].annotate('(NUV-G)$_0$ ='+str(zgthin[0])[:7]+'*T$_{eff}$ '+str(zgthin[1])[:4], xy=(5000, 10))
axes[1, 1].annotate('(NUV-G)$_0$ ='+str(zpthin[0])[:7]+'*T$_{eff}$ '+str(zpthin[1])[:4], xy=(5000, 10))
axes[2, 0].annotate('(NUV-G)$_0$ ='+str(zgthick[0])[:7]+'*T$_{eff}$ '+str(zgthick[1])[:4], xy=(5000, 10))
axes[2, 1].annotate('(NUV-G)$_0$ ='+str(zpthick[0])[:7]+'*T$_{eff}$ '+str(zpthick[1])[:4], xy=(5000, 10))
plt.show()


#######################################################################################
# logg vs tefff plots
#######################################################################################
fig, axes = plt.subplots(2, 2)
axes[0, 0].hist(rcg['LOGG'], bins=20, range=[1.7, 3])
axes[1, 0].hist(rcp['LOGG'], bins=20, range=[1.7, 3])
cmap = axes[0, 1].scatter(rcg['TEFF'], rcg['LOGG'], c=rcg['FE_H'], vmin=-0.5, vmax=0.35, s=1, alpha=0.5)
axes[1, 1].scatter(rcp['TEFF'], rcp['LOGG'], c=rcp['FEH'], vmin=-0.5, vmax=0.35, s=1, alpha=0.5)
axes[0, 1].set_xlim(4200, 5400)
axes[1, 1].set_xlim(4200, 5400)
axes[0, 1].set_ylim(2.1, 2.9)
axes[1, 1].set_ylim(2.1, 2.9)
axes[0, 0].set_xlabel('log g')
axes[1, 0].set_xlabel('log g')
axes[0, 1].set_ylabel('log g')
axes[1, 1].set_ylabel('log g')
axes[0, 1].set_xlabel('T$_{eff}$')
axes[1, 1].set_xlabel('T$_{eff}$')
axes[0, 0].annotate('GAIS', xy=(1.66, 3000))
axes[1, 0].annotate('UVGAPS', xy=(1.66, 109))
cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
plt.show()



fig, axes = plt.subplots(2, 2)
axes[0, 0].hist(r2['LOGG'], bins=20, range=[1, 3])
axes[1, 0].hist(rcg['LOGG'], bins=20, range=[1, 3])
cmap = axes[0, 1].scatter(np.log10(r2['TEFF']), r2['LOGG'], c=r2['FE_H'], vmin=-0.5, vmax=0.35, s=1, alpha=0.1)
axes[1, 1].scatter(np.log10(rcg['TEFF']), rcg['LOGG'], c=rcg['FE_H'], vmin=-0.5, vmax=0.35, s=1, alpha=0.1)
axes[0, 1].set_xlim(3.64, 3.72)
axes[1, 1].set_xlim(3.64, 3.72)
axes[0, 1].set_ylim(2.1, 2.9)
axes[1, 1].set_ylim(2.1, 2.9)
axes[0, 0].set_xlabel('log g')
axes[1, 0].set_xlabel('log g')
axes[0, 1].set_ylabel('log g')
axes[1, 1].set_ylabel('log g')
axes[0, 1].set_xlabel('log T$_{eff}$')
axes[1, 1].set_xlabel('log T$_{eff}$')
axes[0, 0].annotate('strict cut', xy=(1, 4000))
axes[1, 0].annotate('broad cut', xy=(1, 4000))
cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
plt.suptitle('GAIS rc')
plt.show()


#######################################################################################
# feh vs color, plane vs gais
#######################################################################################

x1, y1 = 8.18, -0.28
x2, y2 = 10, 0.57
m = (y2-y1)/(x2-x1)
b = y2 - x2*m
mg = rcg['phot_g_mean_mag']-rcg['distmod']
nuvgrc = rcg['mag_nuv']-rcg['phot_g_mean_mag']
#q = np.where((mg > -0.3) & (nuvgrc > 7) & (nuvgrc < 10) & (m*nuvgrc + b < mg))
#rcg = rcg[q]


xline = np.linspace(6.5, 10, 50)
m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694

thickp,= np.where((rcp['AFE'] > 0.08) & (rcp['AFE'] > (m*rcp['FEH'] + b)))
thinp,= np.where((rcp['AFE'] < 0.08) | (rcp['AFE'] < (m*rcp['FEH'] + b)))
thickg,= np.where((rcg['AFE'] > 0.08) & (rcg['AFE'] > (m*rcg['FEH'] + b)))
thing,= np.where((rcg['AFE'] < 0.08) | (rcg['AFE'] < (m*rcg['FEH'] + b)))
thinp = rcp[thinp]
thickp = rcp[thickp]
thing = rcg[thing]
thickg = rcg[thickg]


# Get fit lines for plot
xg = (rcg['mag_nuv']-rcg['ebv']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebv']*2.85)
xp = (rcp['nuv']-rcp['ebv3d']*7.24)-(rcp['phot_g_mean_mag']-rcp['ebv3d']*2.85)

xgthick = (thickg['mag_nuv']-thickg['ebv']*7.24)-(thickg['phot_g_mean_mag']-thickg['ebv']*2.85)
xpthick = (thickp['nuv']-thickp['ebv3d']*7.24)-(thickp['phot_g_mean_mag']-thickp['ebv3d']*2.85)
xgthin = (thing['mag_nuv']-thing['ebv']*7.24)-(thing['phot_g_mean_mag']-thing['ebv']*2.85)
xpthin = (thinp['nuv']-thinp['ebv3d']*7.24)-(thinp['phot_g_mean_mag']-thinp['ebv3d']*2.85)

yg = rcg['FEH']
yp = rcp['FEH']
ygthick = thickg['FEH']
ypthick = thickp['FEH']
ygthin = thing['FEH']
ypthin = thinp['FEH']

zg, vg = np.polyfit(xg, yg, 1, cov=True)
zp, vp = np.polyfit(xp, yp, 1, cov=True)
pg = np.poly1d(zg)
pp = np.poly1d(zp)

zgthick, vgthick = np.polyfit(xgthick, ygthick, 1, cov=True)
zpthick, vpthick = np.polyfit(xpthick, ypthick, 1, cov=True)
pgthick = np.poly1d(zgthick)
ppthick = np.poly1d(zpthick)

zgthin, vgthin = np.polyfit(xgthin, ygthin, 1, cov=True)
zpthin, vpthin = np.polyfit(xpthin, ypthin, 1, cov=True)
pgthin = np.poly1d(zgthin)
ppthin = np.poly1d(zpthin)

gerr = np.sqrt(np.sum((pg(xg)-yg)**2)/len(xg))
perr = np.sqrt(np.sum((pp(xp)-yp)**2)/len(xp))

gthickerr = np.sqrt(np.sum((pgthick(xgthick)-ygthick)**2)/len(xgthick))
pthickerr = np.sqrt(np.sum((ppthick(xpthick)-ypthick)**2)/len(xpthick))

gthinerr = np.sqrt(np.sum((pgthin(xgthin)-ygthin)**2)/len(xgthin))
pthinerr = np.sqrt(np.sum((ppthin(xpthin)-ypthin)**2)/len(xpthin))



fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
cmap = axes[0, 0].scatter((rcg['mag_nuv']-rcg['ebv']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebv']*2.85), rcg['FEH'], c=rcg['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[0, 1].scatter((rcp['nuv']-rcp['ebv3d']*7.24)-(rcp['phot_g_mean_mag']-rcp['ebv3d']*2.85), rcp['FEH'], c=rcp['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[0, 0].plot(xline, pg(xline), c='red')
axes[0, 1].plot(xline, pp(xline), c='red')

axes[1, 0].scatter((thickg['mag_nuv']-thickg['ebv']*7.24)-(thickg['phot_g_mean_mag']-thickg['ebv']*2.85), thickg['FEH'], c=thickg['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[1, 1].scatter((thickp['nuv']-thickp['ebv3d']*7.24)-(thickp['phot_g_mean_mag']-thickp['ebv3d']*2.85), thickp['FEH'], c=thickp['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[1, 0].plot(xline, pgthick(xline), c='red')
axes[1, 1].plot(xline, pgthick(xline), c='red')

axes[2, 0].scatter((thing['mag_nuv']-thing['ebv']*7.24)-(thing['phot_g_mean_mag']-thing['ebv']*2.85), thing['FEH'], c=thing['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[2, 1].scatter((thinp['nuv']-thinp['ebv3d']*7.24)-(thinp['phot_g_mean_mag']-thinp['ebv3d']*2.85), thinp['FEH'], c=thinp['AFE'], vmin=-0.05, vmax=0.3, s=1, alpha=0.5)
axes[2, 0].plot(xline, pgthin(xline), c='red')
axes[2, 1].plot(xline, pgthin(xline), c='red')


cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('alpha/fe')
axes[0, 0].set_ylabel('[Fe/H]')
axes[1, 0].set_ylabel('[Fe/H]')
axes[2, 0].set_ylabel('[Fe/H]')
axes[2, 0].set_xlabel('(NUV-G)$_0$')
axes[2, 1].set_xlabel('(NUV-G)$_0$')
axes[0, 0].set_xlim(6.5, 10)
axes[0, 1].set_xlim(6.5, 10)
axes[1, 0].set_xlim(6.5, 10)
axes[1, 1].set_xlim(6.5, 10)
axes[2, 0].set_xlim(6.5, 10)
axes[2, 1].set_xlim(6.5, 10)
axes[0, 0].set_ylim(-1.2, 0.6)
axes[0, 1].set_ylim(-1.2, 0.6)
axes[1, 0].set_ylim(-1.2, 0.6)
axes[1, 1].set_ylim(-1.2, 0.6)
axes[2, 0].set_ylim(-1.2, 0.6)
axes[2, 1].set_ylim(-1.2, 0.6)

fig.subplots_adjust(wspace=0, hspace=0)
axes[0, 0].annotate('GAIS', xy=(6.5, 0.4))
axes[0, 1].annotate('UVGAPS', xy=(6.5, 0.4))
axes[1, 0].annotate('Thick', xy=(6.5, 0.4))
axes[2, 0].annotate('Thin', xy=(6.5, 0.4))
axes[0, 0].annotate('[Fe/H] ='+str(zg[0])[:4]+'*(NUV-G)$_0$ '+str(zg[1])[:5], xy=(8.5, -1))
axes[0, 1].annotate('[Fe/H] ='+str(zp[0])[:4]+'*(NUV-G)$_0$ '+str(zp[1])[:5], xy=(8.5, -1))
axes[1, 0].annotate('[Fe/H] ='+str(zgthin[0])[:4]+'*(NUV-G)$_0$ '+str(zgthin[1])[:5], xy=(8.5, -1))
axes[1, 1].annotate('[Fe/H] ='+str(zpthin[0])[:4]+'*(NUV-G)$_0$ '+str(zpthin[1])[:5], xy=(8.5, -1))
axes[2, 0].annotate('[Fe/H] ='+str(zgthick[0])[:4]+'*(NUV-G)$_0$ '+str(zgthick[1])[:5], xy=(8.5, -1))
axes[2, 1].annotate('[Fe/H] ='+str(zpthick[0])[:4]+'*(NUV-G)$_0$ '+str(zpthick[1])[:5], xy=(8.5, -1))
plt.show()




#######################################################################################
# Fe/H vs color, thin and thick disk, gais, by gb
#######################################################################################
xp = np.linspace(6.5, 11.5, 50)
m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
def gb(f, a, b):
	return np.where((f['glat'] > a) & (f['glat'] < b))

'''
# APOGEE
thickg = rcg[np.where((rcg['alphafe'] > 0.08) & (rcg['alphafe'] > (m*rcg['FE_H'] + b)))]
thing = rcg[np.where((rcg['alphafe'] < 0.08) | (rcg['alphafe'] < (m*rcg['FE_H'] + b)))]
brthin = thing['phot_bp_mean_mag']-thing['phot_rp_mean_mag']
ngthin = (thing['mag_nuv']-thing['ebv']*7.24) - (thing['phot_g_mean_mag']-thing['ebv']*2.85)
fehthin = thing['FE_H']
brthick = thickg['phot_bp_mean_mag']-thickg['phot_rp_mean_mag']
ngthick = (thickg['mag_nuv']-thickg['ebv']*7.24) - (thickg['phot_g_mean_mag']-thickg['ebv']*2.85)
fehthick = thickg['FE_H']
'''

# LAMOST
thickg = lg[np.where((lg['AFE'] > 0.08) & (lg['AFE'] > (m*lg['FEH'] + b)))]
thing = lg[np.where((lg['AFE'] < 0.08) | (lg['AFE'] < (m*lg['FEH'] + b)))]
brthin = thing['phot_bp_mean_mag']-thing['phot_rp_mean_mag']
ngthin = (thing['mag_nuv']-thing['ebv']*7.24) - (thing['phot_g_mean_mag']-thing['ebv']*2.85)
fehthin = thing['FEH']
brthick = thickg['phot_bp_mean_mag']-thickg['phot_rp_mean_mag']
ngthick = (thickg['mag_nuv']-thickg['ebv']*7.24) - (thickg['phot_g_mean_mag']-thickg['ebv']*2.85)
fehthick = thickg['FEH']



fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

axes[0, 0].scatter(ngthin[gb(thing,0,10)], fehthin[gb(thing,0,10)], s=4, alpha=0.5, label='Thin')
axes[0, 0].scatter(ngthick[gb(thickg,0,10)], fehthick[gb(thickg,0,10)], s=4, alpha=0.5, label='Thick')
axes[0, 0].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[0, 0].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[0, 1].scatter(ngthin[gb(thing,10, 20)], fehthin[gb(thing,10, 20)], s=4, alpha=0.5)
axes[0, 1].scatter(ngthick[gb(thickg,10, 20)], fehthick[gb(thickg,10, 20)], s=4, alpha=0.5)
axes[0, 1].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[0, 1].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[1, 0].scatter(ngthin[gb(thing,20, 30)], fehthin[gb(thing,20, 30)], s=4, alpha=0.5)
axes[1, 0].scatter(ngthick[gb(thickg,20, 30)], fehthick[gb(thickg,20, 30)], s=4, alpha=0.5)
axes[1, 0].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[1, 0].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[1, 1].scatter(ngthin[gb(thing,30, 40)], fehthin[gb(thing,30, 40)], s=4, alpha=0.5)
axes[1, 1].scatter(ngthick[gb(thickg,30, 40)], fehthick[gb(thickg,30, 40)], s=4, alpha=0.5)
axes[1, 1].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[1, 1].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[2, 0].scatter(ngthin[gb(thing,40, 50)], fehthin[gb(thing,40, 50)], s=4, alpha=0.5)
axes[2, 0].scatter(ngthick[gb(thickg,40, 50)], fehthick[gb(thickg,40, 50)], s=4, alpha=0.5)
axes[2, 0].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[2, 0].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[2, 1].scatter(ngthin[gb(thing,50, 60)], fehthin[gb(thing,50, 60)], s=4, alpha=0.5)
axes[2, 1].scatter(ngthick[gb(thickg,50, 60)], fehthick[gb(thickg,50, 60)], s=4, alpha=0.5)
axes[2, 1].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[2, 1].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[3, 0].scatter(ngthin[gb(thing,60, 70)], fehthin[gb(thing,60, 70)], s=4, alpha=0.5)
axes[3, 0].scatter(ngthick[gb(thickg,60, 70)], fehthick[gb(thickg,60, 70)], s=4, alpha=0.5)
axes[3, 0].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[3, 0].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')

axes[3, 1].scatter(ngthin[gb(thing,70, 90)], fehthin[gb(thing,70, 90)], s=4, alpha=0.5)
axes[3, 1].scatter(ngthick[gb(thickg,70, 90)], fehthick[gb(thickg,70, 90)], s=4, alpha=0.5)
axes[3, 1].plot(xp, 0.21*xp - 1.81, c='green', label='thin eq')
axes[3, 1].plot(xp, 0.27*xp - 2.53, c='red', label='thick eq')


axes[0, 0].legend(loc=4, prop={'size': 10})
axes[3, 0].set_xlabel('(NUV - G)$_0$')
axes[3, 1].set_xlabel('(NUV - G)$_0$')
axes[0, 0].set_ylabel('[Fe/H]')
axes[1, 0].set_ylabel('[Fe/H]')
axes[2, 0].set_ylabel('[Fe/H]')
axes[3, 0].set_ylabel('[Fe/H]')
axes[0, 0].annotate('gb = 0-10', xy=(9., -1.1), size=10)
axes[0, 1].annotate('gb = 10-20', xy=(9.6, -1.1), size=10)
axes[1, 0].annotate('gb = 20-30', xy=(9.6, -1.1), size=10)
axes[1, 1].annotate('gb = 30-40', xy=(9.6, -1.1), size=10)
axes[2, 0].annotate('gb = 40-50', xy=(9.6, -1.1), size=10)
axes[2, 1].annotate('gb = 50-60', xy=(9.6, -1.1), size=10)
axes[3, 0].annotate('gb = 60-70', xy=(9.6, -1.1), size=10)
axes[3, 1].annotate('gb = 70-90', xy=(9.6, -1.1), size=10)
axes[0, 0].set_xlim(6.9, 10.1)
axes[0, 0].set_ylim(-1.2, 0.6)
fig.subplots_adjust(wspace=0, hspace=0)
plt.suptitle('Dust corrected, green=thin, red=thick')





# Now BPRP
fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

axes[0, 0].scatter(brthin[gb(thing,0,10)], fehthin[gb(thing,0,10)], s=4, alpha=0.5, label='Thin')
axes[0, 0].scatter(brthick[gb(thickg,0,10)], fehthick[gb(thickg,0,10)], s=4, alpha=0.5, label='Thick')

axes[0, 1].scatter(brthin[gb(thing,10, 20)], fehthin[gb(thing,10, 20)], s=4, alpha=0.5)
axes[0, 1].scatter(brthick[gb(thicfig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

axes[0, 0].scatter(thing['TEFF'][gb(thing,0,10)], thing['LOGG'][gb(thing,0,10)], s=4, alpha=0.5, label='Thin')
axes[0, 0].scatter(thickg['TEFF'][gb(thickg,0,10)], thickg['LOGG'][gb(thickg,0,10)], s=4, alpha=0.5, label='Thick')

axes[0, 1].scatter(thing['TEFF'][gb(thing,10, 20)], thing['LOGG'][gb(thing,10, 20)], s=4, alpha=0.5)
axes[0, 1].scatter(thickg['TEFF'][gb(thickg,10, 20)], thickg['LOGG'][gb(thickg,10, 20)], s=4, alpha=0.5)

axes[1, 0].scatter(thing['TEFF'][gb(thing,20, 30)], thing['LOGG'][gb(thing,20, 30)], s=4, alpha=0.5)
axes[1, 0].scatter(thickg['TEFF'][gb(thickg,20, 30)], thickg['LOGG'][gb(thickg,20, 30)], s=4, alpha=0.5)

axes[1, 1].scatter(thing['TEFF'][gb(thing,30, 40)], thing['LOGG'][gb(thing,30, 40)], s=4, alpha=0.5)
axes[1, 1].scatter(thickg['TEFF'][gb(thickg,30, 40)], thickg['LOGG'][gb(thickg,30, 40)], s=4, alpha=0.5)

axes[2, 0].scatter(thing['TEFF'][gb(thing,40, 50)], thing['LOGG'][gb(thing,40, 50)], s=4, alpha=0.5)
axes[2, 0].scatter(thickg['TEFF'][gb(thickg,40, 50)], thickg['LOGG'][gb(thickg,40, 50)], s=4, alpha=0.5)

axes[2, 1].scatter(thing['TEFF'][gb(thing,50, 60)], thing['LOGG'][gb(thing,50, 60)], s=4, alpha=0.5)
axes[2, 1].scatter(thickg['TEFF'][gb(thickg,50, 60)], thickg['LOGG'][gb(thickg,50, 60)], s=4, alpha=0.5)

axes[3, 0].scatter(thing['TEFF'][gb(thing,60, 70)], thing['LOGG'][gb(thing,60, 70)], s=4, alpha=0.5)
axes[3, 0].scatter(thickg['TEFF'][gb(thickg,60, 70)], thickg['LOGG'][gb(thickg,60, 70)], s=4, alpha=0.5)

axes[3, 1].scatter(thing['TEFF'][gb(thing,70, 90)], thing['LOGG'][gb(thing,70, 90)], s=4, alpha=0.5)
axes[3, 1].scatter(thickg['TEFF'][gb(thickg,70, 90)], thickg['LOGG'][gb(thickg,70, 90)], s=4, alpha=0.5)


axes[0, 0].legend(loc=4, prop={'size': 10})
axes[3, 0].set_xlabel('TEFF')
axes[3, 1].set_xlabel('TEFF')
axes[0, 0].set_ylabel('log g')
axes[1, 0].set_ylabel('log g')
axes[2, 0].set_ylabel('log g')
axes[3, 0].set_ylabel('log g')
axes[0, 0].annotate('gb = 0-10', xy=(4200, 2.8), size=10)
axes[0, 1].annotate('gb = 10-20', xy=(4200, 2.8), size=10)
axes[1, 0].annotate('gb = 20-30', xy=(4200, 2.8), size=10)
axes[1, 1].annotate('gb = 30-40', xy=(4200, 2.8), size=10)
axes[2, 0].annotate('gb = 40-50', xy=(4200, 2.8), size=10)
axes[2, 1].annotate('gb = 50-60', xy=(4200, 2.8), size=10)
axes[3, 0].annotate('gb = 60-70', xy=(4200, 2.8), size=10)
axes[3, 1].annotate('gb = 70-90', xy=(4200, 2.8), size=10)
axes[0, 0].set_xlim(4200, 5400)
axes[0, 0].set_ylim(2.1, 2.9)
fig.subplots_adjust(wspace=0, hspace=0)
plt.suptitle('GAIS Dust corrected')
kg,10, 20)], fehthick[gb(thickg,10, 20)], s=4, alpha=0.5)

axes[1, 0].scatter(brthin[gb(thing,20, 30)], fehthin[gb(thing,20, 30)], s=4, alpha=0.5)
axes[1, 0].scatter(brthick[gb(thickg,20, 30)], fehthick[gb(thickg,20, 30)], s=4, alpha=0.5)

axes[1, 1].scatter(brthin[gb(thing,30, 40)], fehthin[gb(thing,30, 40)], s=4, alpha=0.5)
axes[1, 1].scatter(brthick[gb(thickg,30, 40)], fehthick[gb(thickg,30, 40)], s=4, alpha=0.5)

axes[2, 0].scatter(brthin[gb(thing,40, 50)], fehthin[gb(thing,40, 50)], s=4, alpha=0.5)
axes[2, 0].scatter(brthick[gb(thickg,40, 50)], fehthick[gb(thickg,40, 50)], s=4, alpha=0.5)

axes[2, 1].scatter(brthin[gb(thing,50, 60)], fehthin[gb(thing,50, 60)], s=4, alpha=0.5)
axes[2, 1].scatter(brthick[gb(thickg,50, 60)], fehthick[gb(thickg,50, 60)], s=4, alpha=0.5)

axes[3, 0].scatter(brthin[gb(thing,60, 70)], fehthin[gb(thing,60, 70)], s=4, alpha=0.5)
axes[3, 0].scatter(brthick[gb(thickg,60, 70)], fehthick[gb(thickg,60, 70)], s=4, alpha=0.5)

axes[3, 1].scatter(brthin[gb(thing,70, 90)], fehthin[gb(thing,70, 90)], s=4, alpha=0.5)
axes[3, 1].scatter(brthick[gb(thickg,70, 90)], fehthick[gb(thickg,70, 90)], s=4, alpha=0.5)


axes[0, 0].legend(loc=4, prop={'size': 10})
axes[3, 0].set_xlabel('BP - RP')
axes[3, 1].set_xlabel('BP - RP')
axes[0, 0].set_ylabel('[Fe/H]')
axes[1, 0].set_ylabel('[Fe/H]')
axes[2, 0].set_ylabel('[Fe/H]')
axes[3, 0].set_ylabel('[Fe/H]')
axes[0, 0].annotate('gb = 0-10', xy=(1.38, -1.1), size=10)
axes[0, 1].annotate('gb = 10-20', xy=(1.45, -1.1), size=10)
axes[1, 0].annotate('gb = 20-30', xy=(1.45, -1.1), size=10)
axes[1, 1].annotate('gb = 30-40', xy=(1.45, -1.1), size=10)
axes[2, 0].annotate('gb = 40-50', xy=(1.45, -1.1), size=10)
axes[2, 1].annotate('gb = 50-60', xy=(1.45, -1.1), size=10)
axes[3, 0].annotate('gb = 60-70', xy=(1.45, -1.1), size=10)
axes[3, 1].annotate('gb = 70-90', xy=(1.45, -1.1), size=10)
axes[0, 0].set_xlim(1.05, 1.5)
axes[0, 0].set_ylim(-1.2, 0.6)
fig.subplots_adjust(wspace=0, hspace=0)




#######################################################################################
# logg vs teff vs gb, thin and thick disk, gais, plane
#######################################################################################
fig, axes = plt.subplots(4, 2, sharex=True, sharey=True)

axes[0, 0].scatter(thing['TEFF'][gb(thing,0,10)], thing['LOGG'][gb(thing,0,10)], s=4, alpha=0.5, label='Thin')
axes[0, 0].scatter(thickg['TEFF'][gb(thickg,0,10)], thickg['LOGG'][gb(thickg,0,10)], s=4, alpha=0.5, label='Thick')

axes[0, 1].scatter(thing['TEFF'][gb(thing,10, 20)], thing['LOGG'][gb(thing,10, 20)], s=4, alpha=0.5)
axes[0, 1].scatter(thickg['TEFF'][gb(thickg,10, 20)], thickg['LOGG'][gb(thickg,10, 20)], s=4, alpha=0.5)

axes[1, 0].scatter(thing['TEFF'][gb(thing,20, 30)], thing['LOGG'][gb(thing,20, 30)], s=4, alpha=0.5)
axes[1, 0].scatter(thickg['TEFF'][gb(thickg,20, 30)], thickg['LOGG'][gb(thickg,20, 30)], s=4, alpha=0.5)

axes[1, 1].scatter(thing['TEFF'][gb(thing,30, 40)], thing['LOGG'][gb(thing,30, 40)], s=4, alpha=0.5)
axes[1, 1].scatter(thickg['TEFF'][gb(thickg,30, 40)], thickg['LOGG'][gb(thickg,30, 40)], s=4, alpha=0.5)

axes[2, 0].scatter(thing['TEFF'][gb(thing,40, 50)], thing['LOGG'][gb(thing,40, 50)], s=4, alpha=0.5)
axes[2, 0].scatter(thickg['TEFF'][gb(thickg,40, 50)], thickg['LOGG'][gb(thickg,40, 50)], s=4, alpha=0.5)

axes[2, 1].scatter(thing['TEFF'][gb(thing,50, 60)], thing['LOGG'][gb(thing,50, 60)], s=4, alpha=0.5)
axes[2, 1].scatter(thickg['TEFF'][gb(thickg,50, 60)], thickg['LOGG'][gb(thickg,50, 60)], s=4, alpha=0.5)

axes[3, 0].scatter(thing['TEFF'][gb(thing,60, 70)], thing['LOGG'][gb(thing,60, 70)], s=4, alpha=0.5)
axes[3, 0].scatter(thickg['TEFF'][gb(thickg,60, 70)], thickg['LOGG'][gb(thickg,60, 70)], s=4, alpha=0.5)

axes[3, 1].scatter(thing['TEFF'][gb(thing,70, 90)], thing['LOGG'][gb(thing,70, 90)], s=4, alpha=0.5)
axes[3, 1].scatter(thickg['TEFF'][gb(thickg,70, 90)], thickg['LOGG'][gb(thickg,70, 90)], s=4, alpha=0.5)


axes[0, 0].legend(loc=4, prop={'size': 10})
axes[3, 0].set_xlabel('TEFF')
axes[3, 1].set_xlabel('TEFF')
axes[0, 0].set_ylabel('log g')
axes[1, 0].set_ylabel('log g')
axes[2, 0].set_ylabel('log g')
axes[3, 0].set_ylabel('log g')
axes[0, 0].annotate('gb = 0-10', xy=(4200, 2.8), size=10)
axes[0, 1].annotate('gb = 10-20', xy=(4200, 2.8), size=10)
axes[1, 0].annotate('gb = 20-30', xy=(4200, 2.8), size=10)
axes[1, 1].annotate('gb = 30-40', xy=(4200, 2.8), size=10)
axes[2, 0].annotate('gb = 40-50', xy=(4200, 2.8), size=10)
axes[2, 1].annotate('gb = 50-60', xy=(4200, 2.8), size=10)
axes[3, 0].annotate('gb = 60-70', xy=(4200, 2.8), size=10)
axes[3, 1].annotate('gb = 70-90', xy=(4200, 2.8), size=10)
axes[0, 0].set_xlim(4200, 5400)
axes[0, 0].set_ylim(2.1, 2.9)
fig.subplots_adjust(wspace=0, hspace=0)
plt.suptitle('GAIS Dust corrected')





rcg = Table.read('rc_gais_gaiadr2_apogee_dust_7-29.txt', format='ascii')      
mg = rcg['phot_g_mean_mag']-rcg['distmod']
nuvgrc = (rcg['mag_nuv']-rcg['ebvsfd']*7.24)-(rcg['phot_g_mean_mag']-rcg['ebvsfd']*2.85)


m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694


q = np.where((mg > -0.3) & (nuvgrc > 7) & (nuvgrc < 10) & (m*nuvgrc + b < mg))
rcg = rcg[q]

l = fits.open('LMRCV1.fits')[1].data
lgal = SkyCoord(l['RA'], l['DEC'], unit=u.deg, frame='icrs')
rcggal = SkyCoord(rcg['RA']*u.deg, rcg['DEC']*u.deg, frame='icrs')
rcgind, lind, angsep, a3 = search_around_sky(rcggal, lgal, 1*u.arcsec)

rcl = hstack([Table(rcg[rcgind]), Table(l[lind])])







m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thickp,= np.where((pl['AFE'] > 0.08) & (pl['AFE'] > (m*pl['FEH'] + b)))
thinp,= np.where((pl['AFE'] < 0.08) | (pl['AFE'] < (m*pl['FEH'] + b)))
thickg,= np.where((lg['AFE'] > 0.08) & (lg['AFE'] > (m*lg['FEH'] + b)))
thing,= np.where((lg['AFE'] < 0.08) | (lg['AFE'] < (m*lg['FEH'] + b)))
#thinp = pl[thinp]
#thickp = pl[thickp]
#thing = gl[thing]
#thickg = gl[thickg]


ngp = (pl['nuv']-pl['ebv3d']*7.24)-(pl['phot_g_mean_mag']-pl['ebv3d']*2.85)
ngg = (lg['mag_nuv']-lg['ebv']*7.24)-(lg['phot_g_mean_mag']-lg['ebv']*2.85)




def mass(f, a, b):
	return np.where((f['MASS'] > a) & (f['MASS'] < b))


fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
cmap = axes[0, 0].scatter(ngg, gl['FEH'], c=gl['MASS'], vmin=0.5, vmax=2, s=1, alpha=0.5)
axes[0, 1].scatter(ngp, pl['FEH'], c=pl['MASS'], vmin=0.5, vmax=2, s=1)
axes[1, 0].scatter(ngg[thickg], gl['FEH'][thickg], c=gl['MASS'][thickg], vmin=0.5, vmax=2, s=1, alpha=0.5)
axes[1, 1].scatter(ngp[thickp], pl['FEH'][thickp], c=pl['MASS'][thickp], vmin=0.5, vmax=2, s=1)
axes[2, 0].scatter(ngg[thing], gl['FEH'][thing], c=gl['MASS'][thing], vmin=0.5, vmax=2, s=1, alpha=0.5)
axes[2, 1].scatter(ngp[thinp], pl['FEH'][thinp], c=pl['MASS'][thinp], vmin=0.5, vmax=2, s=1)

axes[2, 0].set_xlabel('(NUV - G)$_0$')
axes[2, 1].set_xlabel('(NUV - G)$_0$')
axes[0, 0].set_ylabel('[Fe/H]')
axes[1, 0].set_ylabel('[Fe/H]')
axes[2, 0].set_ylabel('[Fe/H]')
axes[0, 0].annotate('GAIS', xy=(9.28, -1.1), size=10)
axes[0, 1].annotate('UVGAPS', xy=(9.7, -1.1), size=10)
axes[1, 0].annotate('Thick', xy=(9.7, -1.1), size=10)
axes[2, 0].annotate('Thin', xy=(9.7, -1.1), size=10)
axes[0, 0].set_xlim(6.5, 10)
axes[0, 0].set_ylim(-1.2, 0.6)
fig.subplots_adjust(wspace=0, hspace=0)

cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('Mass')






fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
cmap = axes[0, 0].scatter(gl['AGE'], gl['AFE'], c=gl['FEH'], vmin=-.8, vmax=.4, s=4)
axes[0, 1].scatter(pl['AGE'], pl['AFE'], c=pl['FEH'], vmin=-.8, vmax=.4, s=4)
axes[1, 0].scatter(gl['AGE'][thickg], gl['AFE'][thickg], c=gl['FEH'][thickg], vmin=-.8, vmax=.4, s=4)
axes[1, 1].scatter(pl['AGE'][thickp], pl['AFE'][thickp], c=pl['FEH'][thickp], vmin=-.8, vmax=.4, s=4)
axes[2, 0].scatter(gl['AGE'][thing], gl['AFE'][thing], c=gl['FEH'][thing], vmin=-.8, vmax=.4, s=4)
axes[2, 1].scatter(pl['AGE'][thinp], pl['AFE'][thinp], c=pl['FEH'][thinp], vmin=-.8, vmax=.4, s=4)



axes[2, 0].set_xlabel('AGE')
axes[2, 1].set_xlabel('AGE')
axes[0, 0].set_ylabel('[Alpha/Fe]')
axes[1, 0].set_ylabel('[Alpha/Fe]')
axes[2, 0].set_ylabel('[Alpha/Fe]')
axes[0, 0].annotate('GAIS', xy=(9.28, -1.1), size=10)
axes[0, 1].annotate('UVGAPS', xy=(9.7, -1.1), size=10)
axes[1, 0].annotate('Thick', xy=(9.7, -1.1), size=10)
axes[2, 0].annotate('Thin', xy=(9.7, -1.1), size=10)
axes[0, 0].set_xlim(0, 14)
axes[0, 0].set_ylim(-0.05, 0.3)
fig.subplots_adjust(wspace=0, hspace=0)

cbar_ax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')














rc = rcg #fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694

rc1 = rc[np.where(rc['FE_H'] < 0.02)]
rc2 = rc[np.where(rc['FE_H'] > 0.02)]


thick1,= np.where((rc1['alphafe'] > 0.08) & (rc1['alphafe'] > (m*rc1['FE_H'] + b)))
thin1,= np.where((rc1['alphafe'] < 0.08) | (rc1['alphafe'] < (m*rc1['FE_H'] + b)))
nuvg1 = ((rc1['mag_nuv']-rc1['ebv3d']*7.24)-(rc1['phot_g_mean_mag']-rc1['ebv3d']*2.85))

thick2,= np.where(rc2['alphafe'] > 0.065)
thin2,= np.where(rc2['alphafe'] < 0.065)
nuvg2 = ((rc2['mag_nuv']-rc2['ebv3d']*7.24)-(rc2['phot_g_mean_mag']-rc2['ebv3d']*2.85))

#rc['alphafe'][afe_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FE_H'])[afe_apo]
#alphafe = rc['alphafe']

# for now use this. err is only for apo
#afeerr = (alphafe * np.sqrt((rc['ALPHA_M_ERR']/rc['ALPHA_M'])**2 + (rc['M_H_ERR']/rc['M_H'])**2 + (rc['FE_H_ERR']/rc['FE_H'])**2))

plt.scatter(rc1['FE_H'][thin1], rc1['alphafe'][thin1], s=20, edgecolor='none', c=nuvg1[thin1], label=r'Low [$\alpha$/Fe]', vmin=7, vmax=11, marker='D', cmap=cm.jet, alpha=0.5)
plt.scatter(rc1['FE_H'][thick1], rc1['alphafe'][thick1], c=nuvg1[thick1], s=20, marker='s', label=r'High [$\alpha$/Fe]', vmin=7, vmax=11, linewidth=1, cmap=cm.jet, alpha=0.5)
plt.scatter(rc2['FE_H'][thin2], rc2['alphafe'][thin2], s=20, edgecolor='none', c=nuvg2[thin2], vmin=7, vmax=11, marker='D', cmap=cm.jet, alpha=0.5)
plt.scatter(rc2['FE_H'][thick2], rc2['alphafe'][thick2], c=nuvg2[thick2], s=20, marker='s', vmin=7, vmax=11, linewidth=1, cmap=cm.jet, alpha=0.5)

xpfeh = np.linspace(-1, 0.02, 50)
plt.plot(xpfeh, xpfeh*m + b, c='black', linewidth=3)
plt.axhline(y=0.065, xmin=0.68, xmax=0.94, c='black', linewidth=3)



plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha$/Fe]')
plt.xlim((-1.2, 0.6))
plt.ylim((-0.1, 0.4))
plt.colorbar().set_label('(NUV - G)$_0$')
leg = plt.legend(scatterpoints=1)
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
plt.show()








'''
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine'))
rc = rc[q]
'''

rc = lg
x1, y1 = 0.02, 0.065
#x2, y2 = -0.88, 0.265
x2, y2 = -0.86, 0.23

#m = (0.265-0.065)/(-0.88-0.02)
#b = 0.0694

m = (y2 - y1) / (x2 - x1)
b = y1 - m*x1

rc1 = rc[np.where(rc['FEH'] < 0.02)]
rc2 = rc[np.where(rc['FEH'] > 0.02)]


thick1,= np.where((rc1['AFE'] > 0.08) & (rc1['AFE'] > (m*rc1['FEH'] + b)))
thin1,= np.where((rc1['AFE'] < 0.08) | (rc1['AFE'] < (m*rc1['FEH'] + b)))
nuvg1 = ((rc1['mag_nuv']-rc1['ebv']*7.24)-(rc1['phot_g_mean_mag']-rc1['ebv']*2.85))

thick2,= np.where(rc2['AFE'] > 0.065)
thin2,= np.where(rc2['AFE'] < 0.065)
nuvg2 = ((rc2['mag_nuv']-rc2['ebv']*7.24)-(rc2['phot_g_mean_mag']-rc2['ebv']*2.85))

#rc['AFE'][AFE_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FEH'])[AFE_apo]
#AFE = rc['AFE']

# for now use this. err is only for apo
#AFEerr = (AFE * np.sqrt((rc['ALPHA_M_ERR']/rc['ALPHA_M'])**2 + (rc['M_H_ERR']/rc['M_H'])**2 + (rc['FEH_ERR']/rc['FEH'])**2))

thin = plt.scatter(rc1['FEH'][thin1], rc1['AFE'][thin1], s=1, edgecolor='none', c=nuvg1[thin1], vmin=7, vmax=11, marker='D', cmap=cm.jet)
thick = plt.scatter(rc1['FEH'][thick1], rc1['AFE'][thick1], edgecolor='none', c=nuvg1[thick1], s=1, marker='D', vmin=7, vmax=11, linewidth=1, cmap=cm.jet)
plt.scatter(rc2['FEH'][thin2], rc2['AFE'][thin2], s=1, edgecolor='none', c=nuvg2[thin2], vmin=7, vmax=11, marker='s', cmap=cm.jet)
plt.scatter(rc2['FEH'][thick2], rc2['AFE'][thick2], c=nuvg2[thick2], edgecolor='none', s=1, marker='s', vmin=7, vmax=11, linewidth=1, cmap=cm.jet)
#plt.errorbar(rc['FEH'], rc['AFE'], xerr=rc['FE_H_err'], ecolor='black', fmt='none', marker='none', mew=0, elinewidth=1.3, **{"zorder":0})



xpfeh = np.linspace(-0.6, 0.02, 50)
plt.plot(xpfeh, xpfeh*m + b, c='black', linewidth=3)
plt.axhline(y=0.065, xmin=0.68, xmax=0.94, c='black', linewidth=3)
plt.axhline(y=0.181, xmin=0.15, xmax=0.335, c='black', linewidth=3)


plt.xlabel('[Fe/H]')
plt.ylabel(r'[$\alpha$/Fe]')
plt.xlim((-1.2, 0.6))
plt.ylim((-0.1, 0.4))
plt.colorbar().set_label('(NUV - G)$_0$')
leg = plt.legend((thin, thick), ('Thin disk', 'Thick disk'), scatterpoints=1)



'''
leg.legendHandles[0].set_color('black')
leg.legendHandles[1].set_color('black')
leg.legendHandles[0]._sizes = [150]
leg.legendHandles[1]._sizes = [150]
'''



x1, y1 = 0.02, 0.065
x2, y2 = -0.86, 0.23
m = (y2 - y1) / (x2 - x1)
b = y1 - m*x1

thickg,= np.where((r['alphafe'] > 0.08) & (r['alphafe'] > (m*r['FE_H'] + b)))
thing,= np.where((r['alphafe'] < 0.08) | (r['alphafe'] < (m*r['FE_H'] + b)))
thing = r[thing]
thickg = r[thickg]
xline = np.linspace(3800, 5300, 50)


# Get fit lines for plot
xg = r['TEFF']
xgthick = thickg['TEFF']
xgthin = thing['TEFF']


yg = (r['nuv_mag']-r['ebv']*7.24)-(r['phot_g_mean_mag']-r['ebv']*2.85)
ygthick = (thickg['nuv_mag']-thickg['ebv']*7.24)-(thickg['phot_g_mean_mag']-thickg['ebv']*2.85)
ygthin = (thing['nuv_mag']-thing['ebv']*7.24)-(thing['phot_g_mean_mag']-thing['ebv']*2.85)



zg, vg = np.polyfit(xg, yg, 1, cov=True)
pg = np.poly1d(zg)
zgthick, vgthick = np.polyfit(xgthick, ygthick, 1, cov=True)
pgthick = np.poly1d(zgthick)
zgthin, vgthin = np.polyfit(xgthin, ygthin, 1, cov=True)
pgthin = np.poly1d(zgthin)
gerr = np.sqrt(np.sum((pg(xg)-yg)**2)/len(xg))
gthickerr = np.sqrt(np.sum((pgthick(xgthick)-ygthick)**2)/len(xgthick))
gthinerr = np.sqrt(np.sum((pgthin(xgthin)-ygthin)**2)/len(xgthin))
