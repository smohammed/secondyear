#WD distances
import numpy as np
from astropy.table import Table
from dustquery import *

# Load all tables
vphas = Table.read('wds_vphasonly.txt', format='ascii')
wd = Table.read('wds_gais_vphas_newgrcut.txt', format='ascii')
vpcut = np.where((vphas['gl'] > 200) & (vphas['gl'] < 250))
wdcut = np.where((wd['gl_galex'] > 200) & (wd['gl_galex'] < 250))
vphas = vphas[vpcut]
wd = wd[wdcut]

# Input constants, set Teff and Radius
sigma = 5.6704 * 10**-5             # erg/cm^2/s/K^4
radius = 637.1 * 10**6              # cm
Temp = np.arange(10, 31, 5)*10**3   # K
Lsun = 3.828 * 10**33               # erg/s
#Lstar = 3.0128 * 10**35             # erg/s, bolometric abs mag
Lum = sigma * Temp**4 * radius**2 * 4 * np.pi   # erg/s
Msun_nuv = 10.18 # mag, from CALSPEC at STScI
Msun_u = 6.44
Msun_g = 5.11
Msun_r = 4.65
Msun_i = 4.54
Mabs_nuv = -2.5 * np.log10(Lum/Lsun) + Msun_nuv
Mabs_u = -2.5 * np.log10(Lum/Lsun) + Msun_u
Mabs_g = -2.5 * np.log10(Lum/Lsun) + Msun_g
Mabs_r = -2.5 * np.log10(Lum/Lsun) + Msun_r
Mabs_i = -2.5 * np.log10(Lum/Lsun) + Msun_i
#Mabs = -2.5 * np.log10(Lum/Lstar)
#freq = 3*10**10/(2.267 * 10**-5)  # NUV
#freq = 3*10**10/(3.543 * 10**-5)  # u
#freq = 3*10**10/(4.770 * 10**-5)  # g
#freq = 3*10**10/(6.231 * 10**-5)  # r
#freq = 3*10**10/(7.625 * 10**-5)  # i

# What table should we use?
comptable = wd
#comptable = vphas[:5000]

# Gather dust data for all WDs
if len(comptable) == len(wd):
    gl = list(comptable['gl_galex'])
    gb = list(comptable['gb_galex'])
if len(comptable) == len(vphas[:5000]):
    gl = list(comptable['gl'])
    gb = list(comptable['gb'])
dust = query(gl, gb, coordsys='gal')
DM = dust['distmod']

# Table will be len(comptable) * len(DM) * len(Mabs) long
table = Table(np.empty((0, 20)))

for i in range(len(comptable)):
    for j in range(len(DM)):
        for k in range(len(Temp)):
            EBV = dust['best'][i][j]
            distmod = DM[j]
            tem = Temp[k]
            # Calculate apparent magnitude model
            mmodel_nuv = -2.5*np.log10(sigma * tem**4 * 4 * np.pi * radius**2 / Lsun) + Msun_nuv + distmod + EBV * 7.76/3.1
            mmodel_u = -2.5*np.log10(sigma * tem**4 * 4 * np.pi * radius**2 / Lsun) + Msun_u + distmod + EBV * 4.239/3.1
            mmodel_g = -2.5*np.log10(sigma * tem**4 * 4 * np.pi * radius**2 / Lsun) + Msun_g + distmod + EBV * 3.303/3.1
            mmodel_r = -2.5*np.log10(sigma * tem**4 * 4 * np.pi * radius**2 / Lsun) + Msun_r + distmod + EBV * 2.285/3.1
            mmodel_i = -2.5*np.log10(sigma * tem**4 * 4 * np.pi * radius**2 / Lsun) + Msun_i + distmod + EBV * 1.698/3.1
            # Compute chi^2
            chi2nuv = (comptable['nuv_mag'][i] - mmodel_nuv)**2
            chi2u = (comptable['u_AB'][i] - mmodel_u)**2
            chi2g = (comptable['g_AB'][i] - mmodel_g)**2
            chi2r = (comptable['r_AB'][i] - mmodel_r)**2
            chi2i = (comptable['i_AB'][i] - mmodel_i)**2
            chi2 = (chi2nuv + chi2u + chi2g + chi2r + chi2i)/5.
            a = [i, comptable['nuv_mag'][i], comptable['u_AB'][i], comptable['g_AB'][i], comptable['r_AB'][i], comptable['i_AB'][i], distmod, tem, EBV, mmodel_nuv, mmodel_u, mmodel_g, mmodel_r, mmodel_i, chi2nuv, chi2u, chi2g, chi2r, chi2i, chi2]
            table.add_row(a)

table.rename_column('col0', 'wdnum')
table.rename_column('col1', 'nuv_obs')
table.rename_column('col2', 'u_obs')
table.rename_column('col3', 'g_obs')
table.rename_column('col4', 'r_obs')
table.rename_column('col5', 'i_obs')
table.rename_column('col6', 'DM')
table.rename_column('col7', 'Temp')
table.rename_column('col8', 'EBV')
table.rename_column('col9', 'nuv_mod')
table.rename_column('col10', 'u_mod')
table.rename_column('col11', 'g_mod')
table.rename_column('col12', 'r_mod')
table.rename_column('col13', 'i_mod')
table.rename_column('col14', 'chi2_nuv')
table.rename_column('col15', 'chi2_u')
table.rename_column('col16', 'chi2_g')
table.rename_column('col17', 'chi2_r')
table.rename_column('col18', 'chi2_i')
table.rename_column('col19', 'chi2_avg')

ascii.write(table, 'wd_GV_model_1ER_Msuncalc.txt', format='basic')


# Find minimum chi^2 value and only plot those
chiind = []
for i in range(0, len(table), 31*5):
    chiind.append(np.argmin(table['chi2'][i:i+31*5])+i)
table = table[chiind]

#cut = np.where(table['chi2'] < 10)
newtable = table  # [cut]
DM = np.unique(newtable['DM'])
temp = np.unique(newtable['Temp'])

chimap = np.empty((len(DM), len(temp)))*0.
chimap[np.isnan(chimap)] = 0

for i in range(len(DM)):
    for j in range(len(temp)):
        chicut = np.where((DM[i] == newtable['DM']) & (temp[j] == newtable['Temp']))
        chimap[i][j] = chimap[i][j] + sum(newtable['chi2_i'][chicut])
        #chimap[i][j] = len(newtable['chi2'][chicut])

yv, xv = np.meshgrid(temp, DM)
plt.pcolormesh(xv, yv, chimap/len(newtable), vmin=0, vmax=1)
plt.xlabel('Distance Modulus')
plt.ylabel('Temperature [K]')
cm = plt.colorbar()
#cm.set_label('$\chi^2/$n$_{bands}$ (m$_{obs}$ - m$_{mod}$)$^2$')
cm.set_label('$\chi^2$ (i$_{obs}$ - i$_{mod}$)$^2$')
#cm.set_label('n$_{WDs}$')
plt.title('GV WD model, 3D dust map, i, 1 R$_{Earth}$')
plt.xlim((4, 19))
plt.show()


# Find minimum chi^2 value and only plot those
chiind = []
for i in range(0, len(table), 31*5):
    chiind.append(np.argmin(table['chi2'][i:i+31*5])+i)
table = table[chiind]


plt.hist2d(table['DM'], table['Temp'])
plt.xlabel('Distance Modulus')
plt.ylabel('Temperature [K]')
cm = plt.colorbar()
cm.set_label('n$_{WDs}$')
plt.title('GV WD model, 0.5 R$_{Earth}$, min $\chi^2$ vals')
plt.show()


old = 0
if old == 1:
    band = 'u'
    DL1 = np.empty((len(comptable), len(Mabs)))
    DL2 = np.empty((len(comptable), len(Mabs)))

    if len(comptable) == len(wd):
        mag = band+'_AB'

    elif len(comptable) == len(vphas[:5000]):
        mag = band+'_AB'

    for i in range(len(comptable)):
        for j in range(len(Mabs)):
            # Distance from distance modulus for each Temp
            DL1[i][j] = 10**((comptable[mag][i]-Mabs[j])/5. + 1)

            # Distance from flux
            DL2[i][j] = np.sqrt(Lum[j]/(4*np.pi * 10**(comptable[mag][i]/-2.5)*3631.*10**-23 * (freq))) / (3.086*10**18)  # parsecs

    DL1 = Table(DL1)
    DL2 = Table(DL2)

    DL1.rename_column('col0', '10k')
    DL1.rename_column('col1', '15k')
    DL1.rename_column('col2', '20k')
    DL1.rename_column('col3', '25k')
    DL1.rename_column('col4', '30k')
    DL2.rename_column('col0', '10k')
    DL2.rename_column('col1', '15k')
    DL2.rename_column('col2', '20k')
    DL2.rename_column('col3', '25k')
    DL2.rename_column('col4', '30k')

    # Gather dust data for all WDs
    if len(comptable) == len(wd):
        gl = list(comptable['gl_galex'])
        gb = list(comptable['gb_galex'])

    if len(comptable) == len(vphas[:5000]):
        gl = list(comptable['gl'])
        gb = list(comptable['gb'])

    dust = query(gl, gb, coordsys='gal')

    # Assign each distance modulus to a dust distance by finding the closest bin
    dist = 10**(np.array(dust['distmod'])/5.+1)  # parsecs

    inddl1 = np.empty((len(comptable), len(Mabs)))
    inddl2 = np.empty((len(comptable), len(Mabs)))

    for i in range(len(dist)-1):
        for j in range(len(comptable)):
            for k in range(len(Mabs)):
                if (DL1[j][k] > dist[i]) and (DL1[j][k] < dist[i+1]):
                    inddl1[j][k] = i
                if (DL2[j][k] > dist[i]) and (DL2[j][k] < dist[i+1]):
                    inddl2[j][k] = i

    inddl1 = Table(inddl1)
    inddl2 = Table(inddl2)
    inddl1['col0'][np.isnan(inddl1['col0'])] = 0
    inddl1['col1'][np.isnan(inddl1['col1'])] = 0
    inddl1['col2'][np.isnan(inddl1['col2'])] = 0
    inddl1['col3'][np.isnan(inddl1['col3'])] = 0
    inddl1['col4'][np.isnan(inddl1['col4'])] = 0
    inddl2['col0'][np.isnan(inddl2['col0'])] = 0
    inddl2['col1'][np.isnan(inddl2['col1'])] = 0
    inddl2['col2'][np.isnan(inddl2['col2'])] = 0
    inddl2['col3'][np.isnan(inddl2['col3'])] = 0
    inddl2['col4'][np.isnan(inddl2['col4'])] = 0

    # Use indices to find E(B-V)
    avdl1 = Table(np.empty((len(comptable), len(Mabs))))
    avdl2 = Table(np.empty((len(comptable), len(Mabs))))
    avdl1.rename_column('col0', '10k')
    avdl1.rename_column('col1', '15k')
    avdl1.rename_column('col2', '20k')
    avdl1.rename_column('col3', '25k')
    avdl1.rename_column('col4', '30k')
    avdl2.rename_column('col0', '10k')
    avdl2.rename_column('col1', '15k')
    avdl2.rename_column('col2', '20k')
    avdl2.rename_column('col3', '25k')
    avdl2.rename_column('col4', '30k')

    temps = ['col0', 'col1', 'col2', 'col3', 'col4']

    for i in range(len(comptable)):
        for j in range(len(temps)):
            avdl1[i][j] = dust['best'][i][int(inddl1[temps[j]][i])]
            avdl2[i][j] = dust['best'][i][int(inddl2[temps[j]][i])]

    # Convert E(B-V) to AV
    avdl1['10k'] = avdl1['10k']*3.1
    avdl1['15k'] = avdl1['15k']*3.1
    avdl1['20k'] = avdl1['20k']*3.1
    avdl1['25k'] = avdl1['25k']*3.1
    avdl1['30k'] = avdl1['30k']*3.1
    avdl2['10k'] = avdl2['10k']*3.1
    avdl2['15k'] = avdl2['15k']*3.1
    avdl2['20k'] = avdl2['20k']*3.1
    avdl2['25k'] = avdl2['25k']*3.1
    avdl2['30k'] = avdl2['30k']*3.1

    # Include which surveys used, the method, what radius used and the band
    if len(comptable) == len(wd):
        ascii.write(DL2, 'wd_GV_dist2_flux_1ER'+band+'.txt', format='basic')
        ascii.write(avdl1, 'wd_GV_av1_dm_1ER_'+band+'.txt', format='basic')
        ascii.write(avdl2, 'wd_GV_av2_flux_1ER_'+band+'.txt', format='basic')

    if len(comptable) == len(vphas[:5000]):
        ascii.write(DL2, 'wd_V_dist2_flux_1ER_'+band+'.txt', format='basic')
        ascii.write(avdl1, 'wd_V_av1_dm_1ER_'+band+'.txt', format='basic')
        ascii.write(avdl2, 'wd_V_av2_flux_1ER_'+band+'.txt', format='basic')

    print band
