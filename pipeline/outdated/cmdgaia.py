from astropy.io import fits
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
from astroML.plotting import scatter_contour
import matplotlib, matplotlib.cm as cm
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 17

#######################################################
# HR diagram by dist and 45 deg slices
#######################################################
ext = 1
sextractor = 1
gais = 0
plane = 0
fuv = 0

if sextractor == 1:
    sgcat = fits.open('../sex_gaia_dust_interp.fits')[1].data
    sggal = SkyCoord(sgcat['gl_sex']*u.deg, sgcat['gb_sex']*u.deg, frame='galactic')

if gais == 1:
    sgcat = fits.open('../gais_gaia_dust.fits')[1].data
    # In the plane
    if plane == 1:
        sgcat = sgcat[np.where((sgcat['gb_gais'] > -10) & (sgcat['gb_gais'] < 10))]

    # Out of the plane
    if plane == 0:
        sgcat = sgcat[np.where((sgcat['gb_gais'] < -10) | (sgcat['gb_gais'] > 10))]

    if fuv == 1:
        sgcat = sgcat[np.where((sgcat['fuv_mag'] != -999.) & (sgcat['fuv_mag'] != -99.) & (sgcat['fuv_mag'] != 99.))]


apo = fits.open('../APOKASKRC_TGAS.fits')[1].data
apogal = SkyCoord(apo['l']*u.deg, apo['b']*u.deg, frame='galactic')
bov = fits.open('../BovyRC_TGAS.fits')[1].data
bovgal = SkyCoord(bov['l']*u.deg, bov['b']*u.deg, frame='galactic')
sg1ind, apoind, angsep1, ang3d = search_around_sky(sggal, apogal, 3*u.arcsec)
sg2ind, bovind, angsep2, ang3d = search_around_sky(sggal, bovgal, 3*u.arcsec)
c1cat = hstack([Table(sgcat)[sg1ind], Table(apo)[apoind]])
c2cat = hstack([Table(sgcat)[sg2ind], Table(bov)[bovind]])
c1cat['angsep'] = angsep1
c2cat['angsep'] = angsep2

sga = sgcat[np.where((sgcat['dist'] > 0) & (sgcat['dist'] < 100) & (sgcat['ebv'] > 0))]
sgb = sgcat[np.where((sgcat['dist'] > 100) & (sgcat['dist'] < 300) & (sgcat['ebv'] > 0))]
sgc = sgcat[np.where((sgcat['dist'] > 300) & (sgcat['dist'] < 600) & (sgcat['ebv'] > 0))]
sgd = sgcat[np.where((sgcat['dist'] > 600) & (sgcat['dist'] < 1000) & (sgcat['ebv'] > 0))]
sge = sgcat[np.where((sgcat['dist'] > 1000) & (sgcat['dist'] < 3000) & (sgcat['ebv'] > 0))]
sgf = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['ebv'] > 0))]

c1a = c1cat[np.where((c1cat['dist'] > 0) & (c1cat['dist'] < 100) & (c1cat['ebv'] > 0))]
c1b = c1cat[np.where((c1cat['dist'] > 100) & (c1cat['dist'] < 300) & (c1cat['ebv'] > 0))]
c1c = c1cat[np.where((c1cat['dist'] > 300) & (c1cat['dist'] < 600) & (c1cat['ebv'] > 0))]
c1d = c1cat[np.where((c1cat['dist'] > 600) & (c1cat['dist'] < 1000) & (c1cat['ebv'] > 0))]
c1e = c1cat[np.where((c1cat['dist'] > 1000) & (c1cat['dist'] < 3000) & (c1cat['ebv'] > 0))]
c1f = c1cat[np.where((c1cat['dist'] > 3000) & (c1cat['ebv'] > 0))]
c2a = c2cat[np.where((c2cat['dist'] > 0) & (c2cat['dist'] < 100) & (c2cat['ebv'] > 0))]
c2b = c2cat[np.where((c2cat['dist'] > 100) & (c2cat['dist'] < 300) & (c2cat['ebv'] > 0))]
c2c = c2cat[np.where((c2cat['dist'] > 300) & (c2cat['dist'] < 600) & (c2cat['ebv'] > 0))]
c2d = c2cat[np.where((c2cat['dist'] > 600) & (c2cat['dist'] < 1000) & (c2cat['ebv'] > 0))]
c2e = c2cat[np.where((c2cat['dist'] > 1000) & (c2cat['dist'] < 3000) & (c2cat['ebv'] > 0))]
c2f = c2cat[np.where((c2cat['dist'] > 3000) & (c2cat['ebv'] > 0))]

cmd = Table.read('../cmdfiles/cmd_merged_zt.txt', format='ascii')

sg = [sga, sgb, sgc, sgd, sge]
c1 = [c1a, c1b, c1c, c1d, c1e]
c2 = [c2a, c2b, c2c, c2d, c2e]

distrange = ['0 < D < 100 pc', '100 < D < 300 pc', '300 < D < 600 pc', '600 < D < 1000 pc', 'D > 1000 pc']
distfilename = ['0-100', '100-300', '300-600', '600-1000', '1000']

horlinelim = [6.5, 6, 4.5, 3.5, 2.5]
verlinelim = [5, 4, 3.3, 3, 3]

zrange = np.unique(cmd['Z'])
agerange = np.unique(cmd['logage'])

colors = ['darkblue', 'blue', 'lightblue', 'green', 'yellow', 'orange', 'red']

for metal in zrange:
    cmd1 = cmd[np.where(cmd['Z'] == metal)]

    for i in range(len(sg)):
        # By gl
        if sextractor == 1:
            sg1 = sg[i][np.where((sg[i]['gl_sex'] > 0) & (sg[i]['gl_sex'] < 45))]
            sg2 = sg[i][np.where((sg[i]['gl_sex'] > 45) & (sg[i]['gl_sex'] < 90))]
            sg3 = sg[i][np.where((sg[i]['gl_sex'] > 90) & (sg[i]['gl_sex'] < 135))]
            sg4 = sg[i][np.where((sg[i]['gl_sex'] > 135) & (sg[i]['gl_sex'] < 180))]
            sg5 = sg[i][np.where((sg[i]['gl_sex'] > 180) & (sg[i]['gl_sex'] < 225))]
            sg6 = sg[i][np.where((sg[i]['gl_sex'] > 225) & (sg[i]['gl_sex'] < 270))]
            sg7 = sg[i][np.where((sg[i]['gl_sex'] > 270) & (sg[i]['gl_sex'] < 315))]
            sg8 = sg[i][np.where((sg[i]['gl_sex'] > 315) & (sg[i]['gl_sex'] < 360))]

        if gais == 1:
            sg1 = sg[i][np.where((sg[i]['gl_gais'] > 0) & (sg[i]['gl_gais'] < 45))]
            sg2 = sg[i][np.where((sg[i]['gl_gais'] > 45) & (sg[i]['gl_gais'] < 90))]
            sg3 = sg[i][np.where((sg[i]['gl_gais'] > 90) & (sg[i]['gl_gais'] < 135))]
            sg4 = sg[i][np.where((sg[i]['gl_gais'] > 135) & (sg[i]['gl_gais'] < 180))]
            sg5 = sg[i][np.where((sg[i]['gl_gais'] > 180) & (sg[i]['gl_gais'] < 225))]
            sg6 = sg[i][np.where((sg[i]['gl_gais'] > 225) & (sg[i]['gl_gais'] < 270))]
            sg7 = sg[i][np.where((sg[i]['gl_gais'] > 270) & (sg[i]['gl_gais'] < 315))]
            sg8 = sg[i][np.where((sg[i]['gl_gais'] > 315) & (sg[i]['gl_gais'] < 360))]

        c11 = c1[i][np.where((c1[i]['l'] > 0) & (c1[i]['l'] < 45))]
        c12 = c1[i][np.where((c1[i]['l'] > 45) & (c1[i]['l'] < 90))]
        c13 = c1[i][np.where((c1[i]['l'] > 90) & (c1[i]['l'] < 135))]
        c14 = c1[i][np.where((c1[i]['l'] > 135) & (c1[i]['l'] < 180))]
        c15 = c1[i][np.where((c1[i]['l'] > 180) & (c1[i]['l'] < 225))]
        c16 = c1[i][np.where((c1[i]['l'] > 225) & (c1[i]['l'] < 270))]
        c17 = c1[i][np.where((c1[i]['l'] > 270) & (c1[i]['l'] < 315))]
        c18 = c1[i][np.where((c1[i]['l'] > 315) & (c1[i]['l'] < 360))]
        c21 = c2[i][np.where((c2[i]['l'] > 0) & (c2[i]['l'] < 45))]
        c22 = c2[i][np.where((c2[i]['l'] > 45) & (c2[i]['l'] < 90))]
        c23 = c2[i][np.where((c2[i]['l'] > 90) & (c2[i]['l'] < 135))]
        c24 = c2[i][np.where((c2[i]['l'] > 135) & (c2[i]['l'] < 180))]
        c25 = c2[i][np.where((c2[i]['l'] > 180) & (c2[i]['l'] < 225))]
        c26 = c2[i][np.where((c2[i]['l'] > 225) & (c2[i]['l'] < 270))]
        c27 = c2[i][np.where((c2[i]['l'] > 270) & (c2[i]['l'] < 315))]
        c28 = c2[i][np.where((c2[i]['l'] > 315) & (c2[i]['l'] < 360))]

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True)

        if ext == 1:
            scatter_contour((sg1['nuv_mag']-sg1['ebv']*7.76)-(sg1['phot_g_mean_mag']-sg1['ebv']*3.303), sg1['Mg']-sg1['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
            scatter_contour((sg2['nuv_mag']-sg2['ebv']*7.76)-(sg2['phot_g_mean_mag']-sg2['ebv']*3.303), sg2['Mg']-sg2['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax2)
            scatter_contour((sg3['nuv_mag']-sg3['ebv']*7.76)-(sg3['phot_g_mean_mag']-sg3['ebv']*3.303), sg3['Mg']-sg3['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax3)
            scatter_contour((sg4['nuv_mag']-sg4['ebv']*7.76)-(sg4['phot_g_mean_mag']-sg4['ebv']*3.303), sg4['Mg']-sg4['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
            scatter_contour((sg5['nuv_mag']-sg5['ebv']*7.76)-(sg5['phot_g_mean_mag']-sg5['ebv']*3.303), sg5['Mg']-sg5['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax5)
            scatter_contour((sg6['nuv_mag']-sg6['ebv']*7.76)-(sg6['phot_g_mean_mag']-sg6['ebv']*3.303), sg6['Mg']-sg6['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax6)
            scatter_contour((sg7['nuv_mag']-sg7['ebv']*7.76)-(sg7['phot_g_mean_mag']-sg7['ebv']*3.303), sg7['Mg']-sg7['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax7)
            scatter_contour((sg8['nuv_mag']-sg8['ebv']*7.76)-(sg8['phot_g_mean_mag']-sg8['ebv']*3.303), sg8['Mg']-sg8['ebv']*3.303, threshold=50, log_counts=True, histogram2d_args=dict(bins=10), plot_args=dict(color='k', markersize=1), contour_args=dict(cmap=cm.gray), ax=ax8)


            ax1.scatter((c11['nuv_mag']-c11['ebv']*7.76)-(c11['phot_g_mean_mag_1']-c11['ebv']*3.303), c11['Mg']-c11['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax2.scatter((c12['nuv_mag']-c12['ebv']*7.76)-(c12['phot_g_mean_mag_1']-c12['ebv']*3.303), c12['Mg']-c12['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax3.scatter((c13['nuv_mag']-c13['ebv']*7.76)-(c13['phot_g_mean_mag_1']-c13['ebv']*3.303), c13['Mg']-c13['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax4.scatter((c14['nuv_mag']-c14['ebv']*7.76)-(c14['phot_g_mean_mag_1']-c14['ebv']*3.303), c14['Mg']-c14['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax5.scatter((c15['nuv_mag']-c15['ebv']*7.76)-(c15['phot_g_mean_mag_1']-c15['ebv']*3.303), c15['Mg']-c15['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax6.scatter((c16['nuv_mag']-c16['ebv']*7.76)-(c16['phot_g_mean_mag_1']-c16['ebv']*3.303), c16['Mg']-c16['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax7.scatter((c17['nuv_mag']-c17['ebv']*7.76)-(c17['phot_g_mean_mag_1']-c17['ebv']*3.303), c17['Mg']-c17['ebv']*3.303, edgecolor='none', c='red', s=80)
            ax8.scatter((c18['nuv_mag']-c18['ebv']*7.76)-(c18['phot_g_mean_mag_1']-c18['ebv']*3.303), c18['Mg']-c18['ebv']*3.303, edgecolor='none', c='red', s=80, label='APO')

            ax1.scatter((c21['nuv_mag']-c21['ebv']*7.76)-(c21['phot_g_mean_mag_1']-c21['ebv']*3.303), c21['Mg']-c21['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax2.scatter((c22['nuv_mag']-c22['ebv']*7.76)-(c22['phot_g_mean_mag_1']-c22['ebv']*3.303), c22['Mg']-c22['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax3.scatter((c23['nuv_mag']-c23['ebv']*7.76)-(c23['phot_g_mean_mag_1']-c23['ebv']*3.303), c23['Mg']-c23['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax4.scatter((c24['nuv_mag']-c24['ebv']*7.76)-(c24['phot_g_mean_mag_1']-c24['ebv']*3.303), c24['Mg']-c24['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax5.scatter((c25['nuv_mag']-c25['ebv']*7.76)-(c25['phot_g_mean_mag_1']-c25['ebv']*3.303), c25['Mg']-c25['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax6.scatter((c26['nuv_mag']-c26['ebv']*7.76)-(c26['phot_g_mean_mag_1']-c26['ebv']*3.303), c26['Mg']-c26['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax7.scatter((c27['nuv_mag']-c27['ebv']*7.76)-(c27['phot_g_mean_mag_1']-c27['ebv']*3.303), c27['Mg']-c27['ebv']*3.303, edgecolor='none', c='blue', s=80)
            ax8.scatter((c28['nuv_mag']-c28['ebv']*7.76)-(c28['phot_g_mean_mag_1']-c28['ebv']*3.303), c28['Mg']-c28['ebv']*3.303, edgecolor='none', c='blue', s=80, label='Bovy')

            ax7.scatter(-100, 100, c='red', label='APO')
            ax7.scatter(-100, 100, c='blue', label='Bovy')
            ax7.legend(scatterpoints=1, loc=2, fontsize='small')


            if sextractor == 1:
                plt.suptitle('SEx+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10')
            if gais == 1:
                if plane == 1:
                    plt.suptitle('GAIS+G, abs(gb) < 10, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10')
                if plane == 0:
                    plt.suptitle('GAIS+G, abs(gb) > 10, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10')

        if ext == 0:
            scatter_contour(sg1['nuv_mag']-sg1['phot_g_mean_mag'], sg1['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax1)
            scatter_contour(sg2['nuv_mag']-sg2['phot_g_mean_mag'], sg2['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax2)
            scatter_contour(sg3['nuv_mag']-sg3['phot_g_mean_mag'], sg3['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax3)
            scatter_contour(sg4['nuv_mag']-sg4['phot_g_mean_mag'], sg4['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax4)
            scatter_contour(sg5['nuv_mag']-sg5['phot_g_mean_mag'], sg5['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax5)
            scatter_contour(sg6['nuv_mag']-sg6['phot_g_mean_mag'], sg6['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax6)
            scatter_contour(sg7['nuv_mag']-sg7['phot_g_mean_mag'], sg7['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax7)
            scatter_contour(sg8['nuv_mag']-sg8['phot_g_mean_mag'], sg8['Mg'], threshold=1000, log_counts=True, histogram2d_args=dict(bins=40), plot_args=dict(color='k', markersize=3), contour_args=dict(cmap=cm.gray), ax=ax8)

            if sextractor == 1:
                plt.suptitle('SEx+G, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10, no extinction')
            if gais == 1:
                if plane == 1:
                    plt.suptitle('GAIS+G, abs(gb) < 10, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10, no extinction')
                if plane == 0:
                    plt.suptitle('GAIS+G, abs(gb) > 10, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10, no extinction')


        for age in range(len(agerange)):
            cmd2 = cmd1[np.where(cmd1['logage'] == agerange[age])]
            ax1.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax2.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax3.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax4.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax5.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax6.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            #ax7.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')
            ax8.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age], linestyle='--')

        if fuv == 0:
            ax1.set_xlim((-1, 11))
            ax1.set_ylim((8, -6))
            ax2.set_xlim((-1, 11))
            ax2.set_ylim((8, -6))
            ax3.set_xlim((-1, 11))
            ax3.set_ylim((8, -6))
            ax4.set_xlim((-1, 11))
            ax4.set_ylim((8, -6))
            ax5.set_xlim((-1, 11))
            ax5.set_ylim((8, -6))
            ax6.set_xlim((-1, 11))
            ax6.set_ylim((8, -6))
            ax7.set_xlim((-1, 11))
            ax7.set_ylim((8, -6))
            ax8.set_xlim((-1, 11))
            ax8.set_ylim((8, -6))
            ax1.annotate('gl = 0-45', xy=(-0.9, 7))
            ax2.annotate('gl = 45-90', xy=(-0.9, 7))
            ax3.annotate('gl = 90-135', xy=(-0.9, 7))
            ax4.annotate('gl = 135-180', xy=(-0.9, 7))
            ax5.annotate('gl = 180-225', xy=(-0.9, 7))
            ax6.annotate('gl = 225-270', xy=(-0.9, 7))
            ax7.annotate('gl = 270-315', xy=(-0.9, 7))
            ax8.annotate('gl = 315-360', xy=(-0.9, 7))
            fig.text(0.5, 0.04, '(NUV - E$_{B-V}$ * 7.76) - (G - E$_{B-V}$ * 3.303)', ha='center')
            fig.text(0.04, 0.5, 'MG - E$_{B-V}$ * 3.303', va='center', rotation='vertical')

        if fuv == 1:
            ax1.set_xlim((-1, 21))
            ax1.set_ylim((7, -5))
            ax2.set_xlim((-1, 21))
            ax2.set_ylim((7, -5))
            ax3.set_xlim((-1, 21))
            ax3.set_ylim((7, -5))
            ax4.set_xlim((-1, 21))
            ax4.set_ylim((7, -5))
            ax5.set_xlim((-1, 21))
            ax5.set_ylim((7, -5))
            ax6.set_xlim((-1, 21))
            ax6.set_ylim((7, -5))
            ax7.set_xlim((-1, 21))
            ax7.set_ylim((7, -5))
            ax8.set_xlim((-1, 21))
            ax8.set_ylim((7, -5))
            ax1.annotate('gl = 0-45', xy=(0, 6))
            ax2.annotate('gl = 45-90', xy=(0, 6))
            ax3.annotate('gl = 90-135', xy=(0, 6))
            ax4.annotate('gl = 135-180', xy=(0, 6))
            ax5.annotate('gl = 180-225', xy=(0, 6))
            ax6.annotate('gl = 225-270', xy=(0, 6))
            ax7.annotate('gl = 270-315', xy=(0, 6))
            ax8.annotate('gl = 315-360', xy=(0, 6))
            fig.text(0.5, 0.04, '(FUV - E$_{B-V}$ * 7.76) - (G - E$_{B-V}$ * 3.303)', ha='center')
            fig.text(0.04, 0.5, 'MG - E$_{B-V}$ * 3.303', va='center', rotation='vertical')

        fig.subplots_adjust(hspace=0, wspace=0)
        '''
        if ext == 1:
            if sextractor == 1:
                plt.savefig('../12-02-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'.png')
            if gais == 1:
                if plane == 1:
                    plt.savefig('../12-02-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_inplane.png')
                if plane == 0:
                    plt.savefig('../12-02-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_outplane.png')

        if ext == 0:
            if sextractor == 1:
                plt.savefig('../12-02-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_noext.png')
            if gais == 1:
                if plane == 1:
                    plt.savefig('../12-02-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_noext_inplane.png')
                if plane == 0:
                    plt.savefig('../12-02-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_noext_outplane.png')
        plt.clf()
        '''
        plt.show()
