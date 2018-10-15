from astropy.io import ascii
import numpy as np
from astropy.table import Table

##########################################################################
# Add LePHARE fix
##########################################################################
#to do:
# 1. Figure out range for each band
# 2. Make array for that range and pull data from Pickles
# 3. Compute: (Sum PICKLE_lambda * LePHARE_lambda) / sum(LePHARE_lambda)
# 4. Make table

c_light = 29979245800  # cm/s
# Open LePHARE files
# Data is set up col1 = wavelength, col2 = flux
directory = '../PICKLES/filt/'
nuvfilt = Table.read(directory+'galex/NUV.pb', format='ascii')

ufilt = Table.read(directory+'sdss/up.pb', format='ascii', data_start=1)
gfilt = Table.read(directory+'sdss/gp.pb', format='ascii', data_start=1)
rfilt = Table.read(directory+'sdss/rp.pb', format='ascii', data_start=1)
ifilt = Table.read(directory+'sdss/ip.pb', format='ascii', data_start=1)
zfilt = Table.read(directory+'sdss/zp.pb', format='ascii', data_start=1)

# Get correct indices for star and filter spectra. 
# Using random star file because their wavelength ranges are all the same
waverange = Table.read('../PICKLES/a0i.sed', format='ascii')
nuv_starind, nuv_filtind = [], []

u_starind, u_filtind = [], []
g_starind, g_filtind = [], []
r_starind, r_filtind = [], []
i_starind, i_filtind = [], []
z_starind, z_filtind = [], []
inda, indn = [], []

for ii in range(len(waverange)):
    for iii in range(len(nuvfilt)):
        if waverange['col1'][ii] == nuvfilt['col1'][iii]:
            nuv_starind.append(ii)
            nuv_filtind.append(iii)

    for iii in range(len(ufilt)):
        if waverange['col1'][ii] == ufilt['col1'][iii]:
            u_starind.append(ii)
            u_filtind.append(iii)
    for iii in range(len(gfilt)):
        if waverange['col1'][ii] == gfilt['col1'][iii]:
            g_starind.append(ii)
            g_filtind.append(iii)
    for iii in range(len(rfilt)):
        if waverange['col1'][ii] == rfilt['col1'][iii]:
            r_starind.append(ii)
            r_filtind.append(iii)
    for iii in range(len(ifilt)):
        if waverange['col1'][ii] == ifilt['col1'][iii]:
            i_starind.append(ii)
            i_filtind.append(iii)
    for iii in range(len(zfilt)):
        if waverange['col1'][ii] == zfilt['col1'][iii]:
            z_starind.append(ii)
            z_filtind.append(iii)

# Extract correct spectra at the correct wavelengths for each band
nuvfilt = nuvfilt[nuv_filtind]

ufilt = ufilt[u_filtind]
gfilt = gfilt[g_filtind]
rfilt = rfilt[r_filtind]
ifilt = ifilt[i_filtind]
zfilt = zfilt[z_filtind]


# Compute sums for denominator: wavelength * S(wavelength)_filter with 3631 Jy
nuvfiltsum = np.sum(nuvfilt['col2'] / nuvfilt['col1']) * 3631 * c_light

ufiltsum = np.sum(ufilt['col2'] / ufilt['col1']) * 3631 * c_light
gfiltsum = np.sum(gfilt['col2'] / gfilt['col1']) * 3631 * c_light
rfiltsum = np.sum(rfilt['col2'] / rfilt['col1']) * 3631 * c_light
ifiltsum = np.sum(ifilt['col2'] / ifilt['col1']) * 3631 * c_light
zfiltsum = np.sum(zfilt['col2'] / zfilt['col1']) * 3631 * c_light

name, nuv, u, g, r, ib, z = [], [], [], [], [], [], []
picklefiles = np.loadtxt('../PICKLES/pickles.txt', dtype='string')

# Compute numerator: wavelength * S()
for i in picklefiles:
    table = Table.read('../PICKLES/'+i, format='ascii')
    name.append(i[:-4])
    nuv.append(np.sum(table['col2'][nuv_starind] * nuvfilt['col2'] * nuvfilt['col1']))

    u.append(np.sum(table['col2'][u_starind] * ufilt['col2'] * ufilt['col1']))
    g.append(np.sum(table['col2'][g_starind] * gfilt['col2'] * gfilt['col1']))
    r.append(np.sum(table['col2'][r_starind] * rfilt['col2'] * rfilt['col1']))
    ib.append(np.sum(table['col2'][i_starind] * ifilt['col2'] * ifilt['col1']))
    z.append(np.sum(table['col2'][z_starind] * zfilt['col2'] * zfilt['col1']))


# Now compute f_lambda
nuv = [x/nuvfiltsum for x in nuv]

u = [x/ufiltsum for x in u]
g = [x/gfiltsum for x in g]
r = [x/rfiltsum for x in r]
ib = [x/ifiltsum for x in ib]
z = [x/zfiltsum for x in z]

pbands = Table([nuv, u, g, r, ib, z], names=('nuv', 'u', 'g', 'r', 'i', 'z'))

'''
# Correction for 0 values
nuvfilt['col2'] = nuvfilt['col2']+0.0000001
Bfilt['col2'] = Bfilt['col2']+0.0000001
Vfilt['col2'] = Vfilt['col2']+0.0000001
Ufilt['col2'] = Ufilt['col2']+0.0000001
Rfilt['col2'] = Rfilt['col2']+0.0000001
Ifilt['col2'] = Ifilt['col2']+0.0000001

Jfilt['col2'] = Jfilt['col2']+0.0000001
Hfilt['col2'] = Hfilt['col2']+0.0000001
Kfilt['col2'] = Kfilt['col2']+0.0000001

ufilt['col2'] = ufilt['col2']+0.0000001
gfilt['col2'] = gfilt['col2']+0.0000001
rfilt['col2'] = rfilt['col2']+0.0000001
ifilt['col2'] = ifilt['col2']+0.0000001
'''

# Compute lambda, only needed for non AB mag
nuvlamb = np.sqrt(np.sum(nuvfilt['col2'])/(np.sum(nuvfilt['col2']*(nuvfilt['col1']+0.0)**-2)))

ulamb = np.sqrt(np.sum(ufilt['col2'])/(np.sum(ufilt['col2']*(ufilt['col1']+0.0)**-2)))
glamb = np.sqrt(np.sum(gfilt['col2'])/(np.sum(gfilt['col2']*(gfilt['col1']+0.0)**-2)))
rlamb = np.sqrt(np.sum(rfilt['col2'])/(np.sum(rfilt['col2']*(rfilt['col1']+0.0)**-2)))
ilamb = np.sqrt(np.sum(ifilt['col2'])/(np.sum(ifilt['col2']*(ifilt['col1']+0.0)**-2)))
zlamb = np.sqrt(np.sum(zfilt['col2'])/(np.sum(zfilt['col2']*(zfilt['col1']+0.0)**-2)))


##########################################################################
# convert to mag
##########################################################################
'''
# these calculations are incorrect, use conversions below these
nuv = -2.5*np.log10(pbands['nuv'] * 3.34*10**4 * (nuvlamb)**2) + 8.9
b = -2.5*np.log10(pbands['b'] * 3.34*10**4 * (blamb)**2) + 8.9
v = -2.5*np.log10(pbands['v'] * 3.34*10**4 * (vlamb)**2) + 8.9
j = -2.5*np.log10(pbands['j']*3.34*10**4 * (Jlamb)**2) + 7.96
h = -2.5*np.log10(pbands['h']*3.34*10**4 * (Hlamb)**2) + 7.526
k = -2.5*np.log10(pbands['k']*3.34*10**4 * (Klamb)**2) + 7.060
u = -2.5*np.log10(pbands['u']*3.34*10**4 * (ulamb)**2) + 8.9lamb
lamb
lamb
g = -2.5*np.log10(pbands['g']*3.34*10**4 * (glamb)**2) + 8.923
r = -2.5*np.log10(pbands['r']*3.34*10**4 * (rlamb)**2) + 9.130
ib = -2.5*np.log10(pbands['i']*3.34*10**4 * (ilamb)**2) + 9.205
'''

nuv = -2.5*np.log10(pbands['nuv'])

u = -2.5*np.log10(pbands['u'])
g = -2.5*np.log10(pbands['g'])
r = -2.5*np.log10(pbands['r'])
ib = -2.5*np.log10(pbands['i'])
z = -2.5*np.log10(pbands['z'])


data = [name, nuv, u, g, r, ib, z]
ascii.write(data, '../picklemags_laphare_10_11_18.txt', names=['name', 'nuv', 'u', 'g', 'r', 'i', 'z'])

