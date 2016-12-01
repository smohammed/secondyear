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
sextractor = 0
gais = 1

if sextractor == 1:
    sgcat = fits.open('sex_gaia_dust_interp.fits')[1].data
if gais == 1:
    sgcat = fits.open('gais_tgasmatch_dust.fits')[1].data
    sgcat = sgcat[np.where((sgcat['gb_gais'] > -10) & (sgcat['gb_gais'] < 10))]

sga = sgcat[np.where((sgcat['dist'] > 0) & (sgcat['dist'] < 100) & (sgcat['ebv'] > 0))]
sgb = sgcat[np.where((sgcat['dist'] > 100) & (sgcat['dist'] < 300) & (sgcat['ebv'] > 0))]
sgc = sgcat[np.where((sgcat['dist'] > 300) & (sgcat['dist'] < 600) & (sgcat['ebv'] > 0))]
sgd = sgcat[np.where((sgcat['dist'] > 600) & (sgcat['dist'] < 1000) & (sgcat['ebv'] > 0))]
sge = sgcat[np.where((sgcat['dist'] > 1000) & (sgcat['dist'] < 3000) & (sgcat['ebv'] > 0))]
sgf = sgcat[np.where((sgcat['dist'] > 3000) & (sgcat['ebv'] > 0))]
cmd = Table.read('cmdfiles/cmd_merged_zt.txt', format='ascii')

sg = [sga, sgb, sgc, sgd, sge]
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
    
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharex=True, sharey=True)

        if ext == 1:
            scatter_contour((sg1['nuv_mag']-sg1['ebv']*7.76)-(sg1['phot_g_mean_mag']-sg1['ebv']*3.303), sg1['Mg']-sg1['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
            scatter_contour((sg2['nuv_mag']-sg2['ebv']*7.76)-(sg2['phot_g_mean_mag']-sg2['ebv']*3.303), sg2['Mg']-sg2['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax2)
            scatter_contour((sg3['nuv_mag']-sg3['ebv']*7.76)-(sg3['phot_g_mean_mag']-sg3['ebv']*3.303), sg3['Mg']-sg3['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax3)
            scatter_contour((sg4['nuv_mag']-sg4['ebv']*7.76)-(sg4['phot_g_mean_mag']-sg4['ebv']*3.303), sg4['Mg']-sg4['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
            scatter_contour((sg5['nuv_mag']-sg5['ebv']*7.76)-(sg5['phot_g_mean_mag']-sg5['ebv']*3.303), sg5['Mg']-sg5['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax5)
            scatter_contour((sg6['nuv_mag']-sg6['ebv']*7.76)-(sg6['phot_g_mean_mag']-sg6['ebv']*3.303), sg6['Mg']-sg6['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax6)
            scatter_contour((sg7['nuv_mag']-sg7['ebv']*7.76)-(sg7['phot_g_mean_mag']-sg7['ebv']*3.303), sg7['Mg']-sg7['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax7)
            scatter_contour((sg8['nuv_mag']-sg8['ebv']*7.76)-(sg8['phot_g_mean_mag']-sg8['ebv']*3.303), sg8['Mg']-sg8['ebv']*3.303,threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax8)

        if ext == 0:
            scatter_contour(sg1['nuv_mag']-sg1['phot_g_mean_mag'], sg1['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax1)
            scatter_contour(sg2['nuv_mag']-sg2['phot_g_mean_mag'], sg2['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax2)
            scatter_contour(sg3['nuv_mag']-sg3['phot_g_mean_mag'], sg3['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax3)
            scatter_contour(sg4['nuv_mag']-sg4['phot_g_mean_mag'], sg4['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax4)
            scatter_contour(sg5['nuv_mag']-sg5['phot_g_mean_mag'], sg5['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax5)
            scatter_contour(sg6['nuv_mag']-sg6['phot_g_mean_mag'], sg6['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax6)
            scatter_contour(sg7['nuv_mag']-sg7['phot_g_mean_mag'], sg7['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax7)
            scatter_contour(sg8['nuv_mag']-sg8['phot_g_mean_mag'], sg8['Mg'],threshold=1000,log_counts=True,histogram2d_args=dict(bins=40),plot_args=dict(color='k',markersize=1), contour_args=dict(cmap=cm.gray), ax=ax8)


        for age in range(len(agerange)):
            cmd2 = cmd1[np.where(cmd1['logage'] == agerange[age])]
            ax1.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax2.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax3.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax4.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax5.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax6.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax7.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
            ax8.plot(cmd2['NUV']-cmd2['G'], cmd2['G'], c=colors[age])
        
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

        if ext == 1:
            if sextractor == 1:
                plt.suptitle('SEx+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10')
            if gais == 1:
                plt.suptitle('GAIS+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10')                
        if ext == 0:
            if sextractor == 1:
                plt.suptitle('SEx+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10, no extinction')
            if gais == 1:
                plt.suptitle('GAIS+G, Ext from S&F11, Hogg parallax, '+distrange[i]+', Z = '+str(metal)+', 8 < log(age/yr) < 10, no extinction')

        fig.subplots_adjust(hspace=0, wspace=0)
        if ext == 1:
            if sextractor == 1:
                plt.savefig('12-01-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'.png')
            if gais == 1:
                plt.savefig('12-01-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'.png')

        if ext == 0:
            if sextractor == 1:
                plt.savefig('12-01-Mgvsnuvg_sex_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_noext.png')
            if gais == 1:
                plt.savefig('12-01-Mgvsnuvg_gais_glcuts_'+distfilename[i]+'pc_Z'+str(metal)+'_noext.png')
        plt.clf()
        #plt.show()  