################################################################
# Calculate errors, get table for paper
################################################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

ra = rc['RA']
dec = rc['Dec']
nuv = rc['nuv_mag'] - rc['ebv']*7.24
nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)
nuverr = rc['nuv_magerr']
Gerr = np.sqrt(rc['phot_g_mean_flux_error']**2 * (-2.5/(np.log(10)*rc['phot_g_mean_flux']))**2)
nuvgerr = np.sqrt(rc['nuv_magerr']**2-Gerr**2)
BP = rc['phot_bp_mean_mag'] - rc['ebv']*2.85
BPRP = (rc['phot_bp_mean_mag']-rc['ebv']*2.85)-(rc['phot_rp_mean_mag']-rc['ebv']*2.85)
BPerr = np.sqrt(rc['phot_bp_mean_flux_error']**2 * (-2.5/(np.log(10)*rc['phot_bp_mean_flux']))**2)
RPerr = np.sqrt(rc['phot_rp_mean_flux_error']**2 * (-2.5/(np.log(10)*rc['phot_rp_mean_flux']))**2)
BPRPerr = np.sqrt(BPerr**2 + RPerr**2)
ebv = rc['ebv']
dm = rc['distmod']
feh = rc['FE_H']
feherr = rc['FE_H_ERR']
teff = rc['TEFF']
alphafe = rc['ALPHAFE']
#afeerr = rc['ALPHAFE_ERR']
afeerr = (alphafe * np.sqrt((rc['ALPHA_M_ERR']/rc['ALPHA_M'])**2 + (rc['M_H_ERR']/rc['M_H'])**2 + (rc['FE_H_ERR']/rc['FE_H'])**2))
afeerr[np.where(afeerr == np.inf)] = 0
afeerr[np.where(afeerr == -1*np.inf)] = 0

ra = ['{:.4f}'.format(x) for x in ra]
dec = ['{:.4f}'.format(x) for x in dec]
nuv = ['{:.2f}'.format(x) for x in nuv]
nuverr = ['{:.2f}'.format(x) for x in nuverr]
nuvg = ['{:.2f}'.format(x) for x in nuvg]
nuvgerr = ['{:.2f}'.format(x) for x in nuvgerr]
BP = ['{:.2f}'.format(x) for x in BP]
BPRP = ['{:.2f}'.format(x) for x in BPRP]
BPerr = ['{:.2f}'.format(x) for x in BPerr]
BPRPerr = ['{:.2f}'.format(x) for x in BPRPerr]
ebv = ['{:.2f}'.format(x) for x in ebv]
dm = ['{:.2f}'.format(x) for x in dm]
feh = ['{:.2f}'.format(x) for x in feh]
feherr = ['{:.2f}'.format(x) for x in feherr]
teff = ['{:.2f}'.format(x) for x in teff]
alphafe = ['{:.2f}'.format(x) for x in alphafe]
afeerr = ['{:.2f}'.format(x) for x in afeerr]

for line in range(len(rc)):
    nuv[line] = nuv[line]+' $\pm$ '+nuverr[line]
    nuvg[line] = nuvg[line]+' $\pm$ '+nuvgerr[line]
    BP[line] = BP[line]+' $\pm$ '+BPerr[line]
    BPRP[line] = BPRP[line]+' $\pm$ '+BPRPerr[line]
    feh[line] = feh[line]+' $\pm$ '+feherr[line]
    alphafe[line] = alphafe[line]+' $\pm$ '+afeerr[line]

table = Table([ra, dec, nuv, nuvg, BP, BPRP, ebv, dm, feh, teff, alphafe])

t1 = table[:10]
ascii.write(t1, Writer=ascii.Latex)

######################################################################
# Color vs TEFF vs Fe/H
######################################################################
m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thick, = np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin, = np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

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

ax1.annotate('(G$_{BP}$ - G$_{RP}$)$_0$ = -0.0004 * T$_{eff}$ + 3.0', xy=(5078, 0.55), color='black', size=13)
ax1.annotate('$\sigma$ = '+str(round(err1, 3)), xy=(5253, 0.62), color='black', size=13)

ax2.annotate('(NUV - G)$_0$ = -0.005 * T$_{eff}$ + 30.8', xy=(5078, 3.2), color='black', size=13)
ax2.annotate('$\sigma$ = '+str(round(err2, 3)), xy=(5253, 3.5), color='black', size=13)

cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
plt.colorbar(cmap, cax=cbar_ax).set_label('[Fe/H]')
plt.show()

##########################################
# Alpha/Fe vs Fe/H vs NUV - G
##########################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine'))
rc = rc[q]

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694

rc1 = rc[np.where(rc['FE_H'] < 0.02)]
rc2 = rc[np.where(rc['FE_H'] > 0.02)]


thick1,= np.where((rc1['ALPHAFE'] > 0.08) & (rc1['ALPHAFE'] > (m*rc1['FE_H'] + b)))
thin1,= np.where((rc1['ALPHAFE'] < 0.08) | (rc1['ALPHAFE'] < (m*rc1['FE_H'] + b)))
nuvg1 = ((rc1['nuv_mag']-rc1['ebv']*7.24)-(rc1['phot_g_mean_mag']-rc1['ebv']*2.85))

thick2,= np.where(rc2['ALPHAFE'] > 0.065)
thin2,= np.where(rc2['ALPHAFE'] < 0.065)
nuvg2 = ((rc2['nuv_mag']-rc2['ebv']*7.24)-(rc2['phot_g_mean_mag']-rc2['ebv']*2.85))

#rc['ALPHAFE'][afe_apo] = (rc['ALPHA_M'] + rc['M_H'] - rc['FE_H'])[afe_apo]
#alphafe = rc['ALPHAFE']

# for now use this. err is only for apo
#afeerr = (alphafe * np.sqrt((rc['ALPHA_M_ERR']/rc['ALPHA_M'])**2 + (rc['M_H_ERR']/rc['M_H'])**2 + (rc['FE_H_ERR']/rc['FE_H'])**2))

plt.scatter(rc1['FE_H'][thin1], rc1['ALPHAFE'][thin1], s=20, edgecolor='none', c=nuvg1[thin1], label=r'Low [$\alpha$/Fe]', vmin=7, vmax=11, marker='D', cmap=cm.jet)
plt.scatter(rc1['FE_H'][thick1], rc1['ALPHAFE'][thick1], c=nuvg1[thick1], s=20, marker='s', label=r'High [$\alpha$/Fe]', vmin=7, vmax=11, linewidth=1, cmap=cm.jet)
plt.scatter(rc2['FE_H'][thin2], rc2['ALPHAFE'][thin2], s=20, edgecolor='none', c=nuvg2[thin2], vmin=7, vmax=11, marker='D', cmap=cm.jet)
plt.scatter(rc2['FE_H'][thick2], rc2['ALPHAFE'][thick2], c=nuvg2[thick2], s=20, marker='s', vmin=7, vmax=11, linewidth=1, cmap=cm.jet)
plt.errorbar(rc['Fe_H'], rc['ALPHAFE'], xerr=rc['FE_H_err'], ecolor='black', fmt='none', marker='none', mew=0, elinewidth=1.3, **{"zorder":0})

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

#################################################
# rc fit lines for Fe/H vs NUV - G vs Alpha/Fe
#################################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

thickfit = rc[thick]
thinfit = rc[thin]

# Key:
# a is not dust corrected. b is.

xa = rc['nuv_mag']-rc['phot_g_mean_mag']
xb = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)
y = rc['FE_H']

za, va = np.polyfit(xa, y, 1, cov=True)
zb, vb = np.polyfit(xb, y, 1, cov=True)
pa = np.poly1d(za)
pb = np.poly1d(zb)

aerr = np.sqrt(np.sum((pa(xa)-y)**2)/len(xa))
berr = np.sqrt(np.sum((pb(xb)-y)**2)/len(xb))

zathin, vathin = np.polyfit(xa[thin], y[thin], 1, cov=True)
zbthin, vbthin = np.polyfit(xb[thin], y[thin], 1, cov=True)
pathin = np.poly1d(zathin)
pbthin = np.poly1d(zbthin)
athinerr = np.sqrt(np.sum((pathin(xa[thin])-y[thin])**2)/len(xa[thin]))
bthinerr = np.sqrt(np.sum((pbthin(xb[thin])-y[thin])**2)/len(xb[thin]))

zathick, vathick = np.polyfit(xa[thick], y[thick], 1, cov=True)
zbthick, vbthick = np.polyfit(xb[thick], y[thick], 1, cov=True)
pathick = np.poly1d(zathick)
pbthick = np.poly1d(zbthick)
athickerr = np.sqrt(np.sum((pathick(xa[thick])-y[thick])**2)/len(xa[thick]))
bthickerr = np.sqrt(np.sum((pbthick(xb[thick])-y[thick])**2)/len(xb[thick]))

xp = np.linspace(6, 10.5, 50)
nuvgerr = np.sqrt(rc['nuv_magerr']**2+rc['Gerr']**2)


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex=True, sharey=True)

# Plotting with no dust correction first
# No dust, all points
cmap = ax1.scatter(xa[thin], y[thin], c=rc['ALPHAFE'][thin], s=20, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax1.errorbar(xa[thin], y[thin], xerr =nuvgerr[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})


ax1.scatter(xa[thick], y[thick], c=rc['ALPHAFE'][thick], s=20, vmin=-0.05, vmax=0.3, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax1.errorbar(xa[thick], y[thick], xerr =nuvgerr[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
ax1.plot(xp, pa(xp), linewidth=4, c='black', zorder=10) 

# Dust, all points
ax2.scatter(xb[thin], y[thin], c=rc['ALPHAFE'][thin], s=20, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax2.errorbar(xb[thin], y[thin], xerr =nuvgerr[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})

ax2.scatter(xb[thick], y[thick], c=rc['ALPHAFE'][thick], s=20, vmin=-0.05, vmax=0.3, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax2.errorbar(xb[thick], y[thick], xerr =nuvgerr[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})

ax2.plot(xp, pb(xp), linewidth=4, c='black', zorder=10)

# No dust, thin disk
ax3.scatter(xa[thin], y[thin], c=rc['ALPHAFE'][thin], s=20, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax3.errorbar(xa[thin], y[thin], xerr =nuvgerr[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
ax3.plot(xp, pathin(xp), linewidth=4, c='black', zorder=10)

# Dust, thin disk
ax4.scatter(xb[thin], y[thin], c=rc['ALPHAFE'][thin], s=20, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax4.errorbar(xb[thin], y[thin], xerr =nuvgerr[thin], yerr=rc['FE_H_ERR'][thin], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
ax4.plot(xp, pbthin(xp), linewidth=4, c='black', zorder=10)

# No dust, thick disk
ax5.scatter(xa[thick], y[thick], c=rc['ALPHAFE'][thick], s=20, vmin=-0.05, vmax=0.3, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax5.errorbar(xa[thick], y[thick], xerr =nuvgerr[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
ax5.plot(xp, pathick(xp), linewidth=4, c='black', zorder=10)

# Dust, thick disk
ax6.scatter(xb[thick], y[thick], c=rc['ALPHAFE'][thick], s=20, vmin=-0.05, vmax=0.3, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
ax6.errorbar(xb[thick], y[thick], xerr =nuvgerr[thick], yerr=rc['FE_H_ERR'][thick], ecolor='black', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
ax6.plot(xp, pbthick(xp), linewidth=4, c='black', zorder=10)

ax5.set_xlabel('NUV - G', fontsize=14)
ax6.set_xlabel('(NUV - G)$_{0}$', fontsize=14)
ax1.set_ylabel('[Fe/H]', fontsize=14)
ax3.set_ylabel('[Fe/H]', fontsize=14)
ax5.set_ylabel('[Fe/H]', fontsize=14)
ax6.legend(scatterpoints=1, loc="lower right")
ax1.set_xlim(5.19, 11)
ax1.set_ylim(-1.1,0.55)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])

ax1.annotate('[Fe/H] = 0.25 * (NUV-G) - 2.22', xy=(7.7, -1.05), color='black', size=13)
ax1.annotate('$\sigma$ = '+str(round(aerr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax2.annotate('[Fe/H] = 0.26 * (NUV-G)$_0$ - 2.20', xy=(7.7, -1.05), color='black', size=13)
ax2.annotate('$\sigma$ = '+str(round(berr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax3.annotate('[Fe/H] = 0.21 * (NUV-G) - 1.81', xy=(7.7, -1.05), color='black', size=13)
ax3.annotate('$\sigma$ = '+str(round(athinerr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax4.annotate('[Fe/H] = 0.22 * (NUV-G)$_0$ - 1.81', xy=(7.7, -1.05), color='black', size=13)
ax4.annotate('$\sigma$ = '+str(round(bthinerr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax5.annotate('[Fe/H] = 0.27 * (NUV-G) - 2.53', xy=(7.7, -1.05), color='black', size=13)
ax5.annotate('$\sigma$ = '+str(round(athickerr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax6.annotate('[Fe/H] = 0.29 * (NUV-G)$_0$ - 2.63', xy=(7.7, -1.05), color='black', size=13)
ax6.annotate('$\sigma$ = '+str(round(bthickerr, 3)), xy=(9.67, -0.92), color='black', size=13)
ax1.annotate('Full sample', xy=(5.2, 0.4), size=13)
ax3.annotate(r'Low [$\alpha$/Fe]', xy=(5.2, 0.4), size=13)
ax5.annotate(r'High [$\alpha$/Fe]', xy=(5.2, 0.4), size=13)
fig.subplots_adjust(hspace=0, wspace=0)
fig.colorbar(cmap, cax=cbar_ax).set_label(r'[$\alpha$/Fe]', fontsize=14)
plt.show()


######################################################################
# FeH vs NUVG vs Alphafe but include alphafe in stats
######################################################################
from scipy.optimize import curve_fit
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)
alphafe = rc['ALPHAFE']

x = [nuvg, alphafe]
xthin = [nuvg[thin], alphafe[thin]]
xthick = [nuvg[thick], alphafe[thick]]
feh = rc['FE_H']

def f(x, a, b, c):
    return a*x[0] + b*x[1] + c

popt, pcov = curve_fit(f, x, feh)
poptthin, pcovthin = curve_fit(f, xthin, feh[thin])
poptthick, pcovthick = curve_fit(f, xthick, feh[thick])

fehp1 = f(x, popt1[0], popt1[1])
sig1 = np.sqrt(np.sum((f1(x, popt1[0], popt1[1])-feh)**2)/len(rc))

err = np.sqrt(np.sum((-y)**2)/len(rc))
err = np.sum(np.sqrt(np.diag(pcov)))
errthin = np.sum(np.sqrt(np.diag(pcovthin)))
errthick = np.sum(np.sqrt(np.diag(pcovthick)))

xp = np.linspace(6, 10.5, 50)
nuvgerr = np.sqrt(rc['nuv_magerr']**2+rc['Gerr']**2)

fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)

# Thick disk
axes[0].scatter(nuvg[thick], feh[thick], c=alphafe[thick], s=10, vmin=-0.05, vmax=0.3, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[0].errorbar(nuvg[thick], feh[thick], xerr=nuvgerr[thick], yerr=rc['FE_H_ERR'][thick], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[0].plot(xp, pthick(xp), linewidth=4, c='red', zorder=10)

# Thin disk
axes[1].scatter(nuvg[thin], feh[thin], c=alphafe[thin], s=10, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[1].errorbar(nuvg[thin], feh[thin], xerr=nuvgerr[thin], yerr=rc['FE_H_ERR'][thin], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[1].plot(xp, pthin(xp), linewidth=4, c='red', zorder=10)

# All
cmap = axes[2].scatter(nuvg, feh, c=alphafe, s=10, vmin=-0.05, vmax=0.3, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[2].errorbar(nuvg, feh, xerr=nuvgerr, yerr=rc['FE_H_ERR'], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[2].plot(xp, p(xp), linewidth=4, c='red', zorder=10)

axes[0].set_ylabel('[Fe/H]', fontsize=14)
axes[1].set_ylabel('[Fe/H]', fontsize=14)
axes[2].set_xlabel('(NUV - G)$_{0}$', fontsize=14)
axes[2].set_ylabel('[Fe/H]', fontsize=14)

fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])

axes[0].annotate('[Fe/H] = 0.16 * (NUV - G)$_{0}$ - 2.04 * r[$\alpha$/Fe] - 1.14', xy=(4400, -1.05), color='black', size=13)
axes[0].annotate('$\sigma$ = '+str(round(thickerr, 3)), xy=(4400, -0.92), color='black', size=13)

axes[1].annotate('[Fe/H] = 0.16 * (NUV - G)$_{0}$ - 2.71 * r[$\alpha$/Fe] - 1.27', xy=(4400, -1.05), color='black', size=13)
axes[1].annotate('$\sigma$ = '+str(round(thinerr, 3)), xy=(4400, -0.92), color='black', size=13)

axes[2].annotate('[Fe/H] = 0.18 * (NUV - G)$_{0}$ - 1.52 * r[$\alpha$/Fe] - 1.47', xy=(4400, -1.05), color='black', size=13)
axes[2].annotate('$\sigma$ = '+str(round(err, 3)), xy=(4400, -0.92), color='black', size=13)

axes[0].annotate('Thick disk', xy=(5150, 0.4), size=13)
axes[1].annotate('Thin disk', xy=(5150, 0.4), size=13)
axes[2].annotate('Full sample', xy=(5150, 0.4), size=13)

fig.subplots_adjust(hspace=0, wspace=0)
fig.colorbar(cmap, cax=cbar_ax).set_label(r'[$\alpha$/Fe]', fontsize=14)
plt.show()

######################################################################
# FeH as a function of Teff with NUVG colorbar
######################################################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

thickfit = rc[thick]
thinfit = rc[thin]

x = rc['TEFF']
y = rc['FE_H']

z, v = np.polyfit(x, y, 1, cov=True)
p = np.poly1d(z)
err = np.sqrt(np.sum((p(x)-y)**2)/len(rc))

zthin, vthin = np.polyfit(x[thin], y[thin], 1, cov=True)
pthin = np.poly1d(zthin)
thinerr = np.sqrt(np.sum((pthin(x[thin])-y[thin])**2)/len(rc[thin]))

zthick, vthick = np.polyfit(x[thick], y[thick], 1, cov=True)
pthick = np.poly1d(zthick)
thickerr = np.sqrt(np.sum((pthick(x[thick])-y[thick])**2)/len(rc[thick]))

xp = np.linspace(4500, 5200, 50)
nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)


fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)

# Thick disk
axes[2].scatter(x[thick], y[thick], c=nuvg[thick], s=10, vmin=7, vmax=11, marker='s', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[2].errorbar(x[thick], y[thick], xerr=rc['TEFF_ERR'][thick], yerr=rc['FE_H_ERR'][thick], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[2].plot(xp, pthick(xp), linewidth=4, c='red', zorder=10)

# Thin disk
axes[1].scatter(x[thin], y[thin], c=nuvg[thin], s=10, vmin=7, vmax=11, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[1].errorbar(x[thin], y[thin], xerr=rc['TEFF_ERR'][thin], yerr=rc['FE_H_ERR'][thin], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[1].plot(xp, pthin(xp), linewidth=4, c='red', zorder=10)

# All
cmap = axes[0].scatter(x, y, c=nuvg, s=10, vmin=7, vmax=11, marker='D', edgecolor='black', cmap=cm.plasma_r, **{"zorder":5})
axes[0].errorbar(x, y, xerr=rc['TEFF_ERR'], yerr=rc['FE_H_ERR'], ecolor='gray', fmt='None', marker='None', mew=0, elinewidth=1.3, **{"zorder":0})
axes[0].plot(xp, p(xp), linewidth=4, c='red', zorder=10)

axes[0].set_ylabel('[Fe/H]', fontsize=14)
axes[1].set_ylabel('[Fe/H]', fontsize=14)
axes[2].set_xlabel('T$_{eff}$', fontsize=14)
axes[2].set_ylabel('[Fe/H]', fontsize=14)

axes[0].set_xlim(4400, 5300)
axes[0].set_ylim(-1.1, 0.55)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])

axes[2].annotate('[Fe/H] = -1.59 * 10$^{-3}$ * T$_{eff}$ + 7.34', xy=(4400, -1.05), color='black', size=13)
axes[2].annotate('$\sigma$ = '+str(round(thickerr, 3)), xy=(4400, -0.92), color='black', size=13)

axes[1].annotate('[Fe/H] = -1.16 * 10$^{-3}$ * T$_{eff}$ + 5.51', xy=(4400, -1.05), color='black', size=13)
axes[1].annotate('$\sigma$ = '+str(round(thinerr, 3)), xy=(4400, -0.92), color='black', size=13)

axes[0].annotate('[Fe/H] = -1.26 * 10$^{-3}$ * T$_{eff}$ + 5.94', xy=(4400, -1.05), color='black', size=13)
axes[0].annotate('$\sigma$ = '+str(round(err, 3)), xy=(4400, -0.92), color='black', size=13)

axes[2].annotate('Thick disk', xy=(5150, 0.4), size=13)
axes[1].annotate('Thin disk', xy=(5150, 0.4), size=13)
axes[0].annotate('Full sample', xy=(5150, 0.4), size=13)

fig.subplots_adjust(hspace=0, wspace=0)
fig.colorbar(cmap, cax=cbar_ax).set_label('(NUV - G)$_0$', fontsize=14)
plt.show()

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

######################################################################
# RC paper histograms
######################################################################
rc = fits.open('asc_gaia_aporc_match_dust-05-03-18.fits')[1].data
q = np.where((rc['ebv'] > 0) & (rc['Fe_H_err'] > 0) & (rc['phot_bp_mean_mag'] > 0) & (rc['phot_rp_mean_mag'] > 0) & (rc['Classification'] == 'RC_Pristine') & (rc['dist'] < 3500) & (rc['visibility_periods_used'] > 8) & (rc['parallax_error']/rc['parallax'] < 0.1))
rc = rc[q]

m = (0.265-0.065)/(-0.88-0.02)
b = 0.0694
thick,= np.where((rc['ALPHAFE'] > 0.08) & (rc['ALPHAFE'] > (m*rc['FE_H'] + b)))
thin,= np.where((rc['ALPHAFE'] < 0.08) | (rc['ALPHAFE'] < (m*rc['FE_H'] + b)))

thin = rc[thin]
thick = rc[thick]

fig, axes = plt.subplots(2, 2)
axes[0, 0].hist(rc['gb_asc'], range=[-90,90])
axes[0, 0].hist(thin['gb_asc'], range=[-90,90], histtype='step',color='red', fill=False, stacked=True)
axes[0, 0].hist(thick['gb_asc'], range=[-90,90], histtype='step', color='black', fill=False, stacked=True)

axes[0, 1].hist(np.log10(rc['dist']), range=[1.8,3.5])
axes[0, 1].hist(np.log10(thin['dist']), range=[1.8,3.5], histtype='step', color='red', stacked=True, fill=False, label=r'Low [$\alpha$/Fe]')
axes[0, 1].hist(np.log10(thick['dist']), range=[1.8,3.5], histtype='step', color='black', stacked=True, fill=False, label=r'High [$\alpha$/Fe]')

axes[1, 0].hist(rc['nuv_mag'], range=[12.5, 22.5])
axes[1, 0].hist(thin['nuv_mag'],range=[12.5, 22.5], histtype='step', color='red', stacked=True, fill=False)
axes[1, 0].hist(thick['nuv_mag'],range=[12.5, 22.5], histtype='step', color='black', stacked=True, fill=False)

axes[1, 1].hist(np.log10(rc['ebv']), range=[-3,0])
axes[1, 1].hist(np.log10(thin['ebv']), range=[-3,0], histtype='step', color='red', stacked=True, fill=False)
axes[1, 1].hist(np.log10(thick['ebv']), range=[-3,0], histtype='step', color='black', stacked=True, fill=False)

axes[0, 0].set_xlabel('Galactic Latitude')
axes[0, 1].set_xlabel('log Distance [pc]')
axes[1, 0].set_xlabel('NUV')
axes[1, 1].set_xlabel('log E(B-V)')
axes[0, 1].legend()
plt.show()

# Color excess
ex = rc[np.where(nuvg*m + b < rc['Fe_H'])]
m = (-0.4 - 0.2)/(6 - 7)
b = -4
fig, axes = plt.subplots(2, 2)
axes[0, 0].hist(rc['gb'], range=[-90,90])
axes[0, 0].hist(ex['gb'], range=[-90,90], histtype='step',color='red', fill=False, stacked=True)

axes[0, 1].hist(np.log10(rc['dist']), range=[1.8,3.5])
axes[0, 1].hist(np.log10(ex['dist']), range=[1.8,3.5], histtype='step', color='red', stacked=True, fill=False, label='NUVG Excess')

axes[1, 0].hist(rc['nuv'], range=[12.5, 22.5])
axes[1, 0].hist(ex['nuv'],range=[12.5, 22.5], histtype='step', color='red', stacked=True, fill=False)

axes[1, 1].hist(np.log10(rc['ebv']), range=[-3,0])
axes[1, 1].hist(np.log10(ex['ebv']), range=[-3,0], histtype='step', color='red', stacked=True, fill=False)

axes[0, 0].set_xlabel('Galactic Latitude')
axes[0, 1].set_xlabel('log(Distance) [pc]')
axes[1, 0].set_xlabel('NUV')
axes[1, 1].set_xlabel('log(E(B-V))')
axes[0, 1].legend()
plt.show()

######################################################################
# Getting stats for test sample
######################################################################
gg = fits.open('galex-asc-gaia-match_smaller.fits', memmap=True)[1].data
pc = 1000./gg['parallax']
negpar = np.where((pc > 0) & (gg['visibility_periods_used'] > 8) & (gg['parallax_error']/gg['parallax'] < 0.1) & (gg['phot_bp_mean_mag'] > 0) & (gg['phot_rp_mean_mag'] > 0)) 
ra = gg['gaia_ra']
dec = gg['gaia_dec']
nuv = gg['mag_nuv']
G = gg['phot_g_mean_mag']
bp = gg['phot_bp_mean_mag']
rp = gg['phot_rp_mean_mag']

ra = ra[negpar]
dec = dec[negpar]
pc = pc[negpar]
nuv = nuv[negpar]
G = G[negpar]
bp = bp[negpar]
rp = rp[negpar]
distmod = 5. * np.log10(pc) - 5.
MG = G - distmod
nuvg = nuv - G

rccut = np.where((MG < 0.9) & (MG > -0.1) & (nuvg < 11.5) & (nuvg > 6))
ra = ra[rccut]
dec = dec[rccut]
pc = pc[rccut]
nuv = nuv[rccut]
G = G[rccut]
bp = bp[rccut]
rp = rp[rccut]
distmod = distmod[rccut]
MG = G - distmod
nuvg = nuv - G

table = Table([ra, dec, nuv, G, bp, rp, distmod, MG, pc, nuvg], names=['ra', 'dec', 'nuv', 'G', 'bp', 'rp', 'distmod', 'MG', 'pc', 'nuvg'])
ascii.write(table, 'gais_rcbox_testset.txt', format='basic')

from dustmaps.bayestar import BayestarQuery
coords = SkyCoord(ra*u.deg, dec*u.deg, distance=pc*u.pc, frame='icrs')
bayestar = BayestarQuery()
ebv = bayestar(coords, mode='median')
table['ebv'] = ebv
ascii.write(table, 'gais_rcbox_testset.txt', format='basic')

#then combine with apo, then get alphafe
table['ALPHAFE'] = table['ALPHA_M'] + table['M_H'] - table['FE_H']

plt.scatter(nuvg, comb['FE_H'], c=alphafe, s=1, vmin=-0.5, vmax=0.3, label='Test set')
plt.plot(xp, 0.26*xp - 2.2, label='Training set fit')
plt.plot(xp, p[0]*xp + p[1], label='Test set fit')
plt.xlabel('(NUV - G)$_0$')
plt.ylabel('[Fe/H]')
plt.xlim(6,11)
plt.ylim(-1, 0.5)
plt.legend(scatterpoints=1, loc=4)
plt.colorbar().set_label(r'[$\alpha$/Fe]')
plt.show()

#####################################################################
# Fe/H comparisons
#####################################################################
nuvg = (rc['nuv_mag']-rc['ebv']*7.24)-(rc['phot_g_mean_mag']-rc['ebv']*2.85)
alphafe = rc['ALPHAFE']
x = [nuvg, alphafe]
feh = rc['FE_H']

def f1(x, a, b):
    return a*x[0] + b
popt1, pcov1 = curve_fit(f1, x, feh)
fehp1 = f1(x, popt1[0], popt1[1])
sig1 = np.sqrt(np.sum((f1(x, popt1[0], popt1[1])-feh)**2)/len(rc))

def f2(x, a, b, c):
    return a*x[0] + b*x[0]**2 + c
popt2, pcov2 = curve_fit(f2, x, feh)
fehp2 = f2(x, popt2[0], popt2[1], popt2[2])
sig2 = np.sqrt(np.sum((fehp2-feh)**2)/len(rc))

def f3(x, a, b, c):
    return a*x[0] + b*x[1] + c
popt3, pcov3 = curve_fit(f3, x, feh)
fehp3 = f3(x, popt3[0], popt3[1], popt3[2])
sig3 = np.sqrt(np.sum((fehp3-feh)**2)/len(rc))

def f4(x, a, b, c, d):
    return a*x[0] + b*x[1] + c*x[0]**2 + d
popt4, pcov4 = curve_fit(f4, x, feh)
fehp4 = f4(x, popt4[0], popt4[1], popt4[2], popt4[3])
sig4 = np.sqrt(np.sum((fehp4-feh)**2)/len(rc))

def f5(x, a, b, c, d):
    return a*x[0] + b*x[1] + c*x[1]**2 + d
popt5, pcov5 = curve_fit(f5, x, feh)
fehp5 = f5(x, popt5[0], popt5[1], popt5[2], popt5[3])
sig5 = np.sqrt(np.sum((fehp5-feh)**2)/len(rc))

def f6(x, a, b, c, d, e):
    return a*x[0] + b*x[1] + c*x[0]**2 + d*x[1]**2 + e
popt6, pcov6 = curve_fit(f6, x, feh)
fehp6 = f6(x, popt6[0], popt6[1], popt6[2], popt6[3], popt6[4])
sig6 = np.sqrt(np.sum((fehp6-feh)**2)/len(rc))

fig, axes = plt.subplots(3, 2, sharex=True, sharey=True)
cmap = axes[0, 0].scatter(feh, feh-fehp1, s=1, c=nuvg, vmin=7, vmax=11)
axes[0, 0].axhline(y=0, c='black')
axes[0, 0].annotate('[Fe/H]$_{phot}$ = '+str(round(popt1[0], 2))+'*nuvg + '+str(round(popt1[1], 2))+', $\sigma$='+str(round(sig1, 3)), xy=(-1.2, -0.7), color='black', size=10)
axes[0, 0].set_xlim(-1.2, 0.6)
axes[0, 0].set_ylim(-0.75, 0.75)

axes[0, 1].scatter(feh, feh-fehp2, s=1, c=nuvg, vmin=7, vmax=11)
axes[0, 1].axhline(y=0, color='black')
axes[0, 1].annotate('[Fe/H]$_{phot}$ = '+str(round(popt2[0], 2))+'*nuvg + '+str(round(popt2[1], 2))+'*nuvg$^2$ + '+str(round(popt2[2], 2))+', $\sigma$='+str(round(sig2, 3)), xy=(-1.2, -0.7), color='black', size=10)

axes[1, 0].scatter(feh, feh-fehp3, s=1, c=nuvg, vmin=7, vmax=11)
axes[1, 0].axhline(y=0, color='black')
axes[1, 0].annotate('[Fe/H]$_{phot}$ = '+str(round(popt3[0], 2))+'*nuvg + '+str(round(popt3[1], 2))+'*alphafe + '+str(round(popt3[2], 2))+', $\sigma$='+str(round(sig3, 3)), xy=(-1.2, -0.7), color='black', size=10)
axes[1, 0].set_ylabel('$\Delta$ [Fe/H]')

axes[1, 1].scatter(feh, feh-fehp4, s=1, c=nuvg, vmin=7, vmax=11)
axes[1, 1].axhline(y=0, color='black')
axes[1, 1].annotate('[Fe/H]$_{phot}$ = '+str(round(popt4[0], 2))+'*nuvg + '+str(round(popt4[1], 2))+'*alphafe + '+str(round(popt4[2], 2))+'*nuvg$^2$ + '+str(round(popt4[3], 2))+', $\sigma$='+str(round(sig4, 3)), xy=(-1.2, -0.7), color='black', size=10)

axes[2, 0].scatter(feh, feh-fehp5, s=1, c=nuvg, vmin=7, vmax=11)
axes[2, 0].axhline(y=0, color='black')
axes[2, 0].annotate('[Fe/H]$_{phot}$ = '+str(round(popt5[0], 2))+'*nuvg + '+str(round(popt5[1], 2))+'*alphafe + '+str(round(popt5[2], 2))+'*alphafe$^2$ + '+str(round(popt5[3], 2))+', $\sigma$='+str(round(sig5, 3)), xy=(-1.2, -0.7), color='black', size=10)
axes[2, 0].set_xlabel('Spec [Fe/H]')

axes[2, 1].scatter(feh, feh-fehp6, s=1, c=nuvg, vmin=7, vmax=11)
axes[2, 1].axhline(y=0, color='black')
axes[2, 1].annotate('[Fe/H]$_{phot}$ = '+str(round(popt6[0], 2))+'*nuvg + '+str(round(popt6[1], 2))+'*alphafe + '+str(round(popt6[2], 2))+'*nuvg$^2$ + '+str(round(popt6[3], 2))+'*alphafe$^2$ + '+str(round(popt6[4], 2))+', $\sigma$='+str(round(sig6, 3)), xy=(-1.2, -0.7), color='black', size=10)
axes[2, 1].set_xlabel('Spec [Fe/H]')

fig.subplots_adjust(hspace=0, wspace=0)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
fig.colorbar(cmap, cax=cbar_ax).set_label('NUV - G', fontsize=15)
fig.suptitle('Thick only')

plt.show()


#####################################################################
# Fe/H comparison for paper
#####################################################################
def f3(x, a, b, c):
    return a*x[0] + b*x[1] + c
x = [nuvg, afe]
popt3, pcov3 = curve_fit(f3, x, feh)
fehp3 = f3(x, popt3[0], popt3[1], popt3[2])
err3 = np.sqrt(np.sum((fehp3-feh)**2)/len(rc))

def f3thin(x, a, b, c):
    return a*x[0] + b*x[1] + c
xthin = [nuvg[thin], afe[thin]]
popt3thin, pcov3thin = curve_fit(f3thin, xthin, feh[thin])
fehp3thin = f3thin(xthin, popt3thin[0], popt3thin[1], popt3thin[2])
err3thin = np.sqrt(np.sum((fehp3thin-feh[thin])**2)/len(rc[thin]))

def f3thick(x, a, b, c):
    return a*x[0] + b*x[1] + c
xthick = [nuvg[thick], afe[thick]]
popt3thick, pcov3thick = curve_fit(f3thick, xthick, feh[thick])
fehp3thick = f3thick(xthick, popt3thick[0], popt3thick[1], popt3thick[2])
err3thick = np.sqrt(np.sum((fehp3thick-feh[thick])**2)/len(rc[thick]))


fig, axes = plt.subplots(3, 1, sharex=True, sharey=True)
cmap = axes[0].scatter(feh, feh-fehp3, s=1, c=nuvg, vmin=7, vmax=11)
axes[0].axhline(y=0, color='black')
axes[0].annotate('[Fe/H]$_{phot}$ = '+str(round(popt3[0], 2))+r' $\times$ (NUV - G)$_0$ + '+str(round(popt3[1], 2))+r' $\times$ [$\alpha$/Fe] + '+str(round(popt3[2], 2))+', $\sigma$='+str(round(err3, 3)), xy=(-1.2, -0.7), color='black', size=13)
axes[0].set_ylabel('$\Delta$ [Fe/H]')

axes[1].scatter(feh[thin], feh[thin]-fehp3thin, s=1, c=nuvg[thin], vmin=7, vmax=11)
axes[1].axhline(y=0, color='black')
axes[1].annotate('[Fe/H]$_{phot}$ = '+str(round(popt3thin[0], 2))+r' $\times$ (NUV - G)$_0$ + '+str(round(popt3thin[1], 2))+r' $\times$ [$\alpha$/Fe] + '+str(round(popt3thin[2], 2))+', $\sigma$='+str(round(err3thin, 3)) + r', Low [$\alpha$/Fe]', xy=(-1.2, -0.7), color='black', size=13)
axes[1].set_ylabel('$\Delta$ [Fe/H]')


axes[2].scatter(feh[thick], feh[thick]-fehp3thick, s=1, c=nuvg[thick], vmin=7, vmax=11)
axes[2].axhline(y=0, color='black')
axes[2].annotate('[Fe/H]$_{phot}$ = '+str(round(popt3thick[0], 2))+r' $\times$ (NUV - G)$_0$ + '+str(round(popt3thick[1], 2))+r' $\times$ [$\alpha$/Fe] + '+str(round(popt3thick[2], 2))+', $\sigma$='+str(round(err3thick, 3)) + r', High [$\alpha$/Fe]', xy=(-1.2, -0.7), color='black', size=13)
axes[2].set_ylabel('$\Delta$ [Fe/H]')
axes[2].set_xlabel('[Fe/H]')
fig.subplots_adjust(hspace=0, wspace=0)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(right=.84)
cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
fig.colorbar(cmap, cax=cbar_ax).set_label('(NUV - G)$_0$', fontsize=15)
plt.show()



#####################################################################
# log g plots
#####################################################################
fig, axes = plt.subplots(1, 2)
cmap = axes[1].scatter(np.log10(rc['TEFF']), rc['logg'], s=5, c=rc['FE_H'], vmin=-0.5, vmax=0.35)

axes[0].hist(rc['logg'], bins=20)
axes[0].set_xlabel('log g')
axes[1].set_ylabel('log g')
axes[1].set_xlabel('log T$_{eff}$')
plt.colorbar(cmap).set_label('[Fe/H]')
plt.show()
