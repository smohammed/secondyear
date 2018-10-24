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
directory = '../PICKLES/filt/'
nuvfilt = Table.read(directory+'galex/NUV.pb', format='ascii')
Bfilt = Table.read(directory+'wfi/B.pb', format='ascii', data_start=1)
Vfilt = Table.read(directory+'wfi/V.pb', format='ascii', data_start=1)
Ufilt = Table.read(directory+'wfi/U.pb', format='ascii', data_start=1)
Rfilt = Table.read(directory+'wfi/R.pb', format='ascii', data_start=1)
Ifilt = Table.read(directory+'wfi/I.pb', format='ascii', data_start=1)

Jfilt = Table.read(directory+'2mass/J.pb', format='ascii', data_start=1)
Hfilt = Table.read(directory+'2mass/H.pb', format='ascii', data_start=1)
Kfilt = Table.read(directory+'2mass/Ks.pb', format='ascii', data_start=1)

ufilt = Table.read(directory+'sdss/up.pb', format='ascii', data_start=1)
gfilt = Table.read(directory+'sdss/gp.pb', format='ascii', data_start=1)
rfilt = Table.read(directory+'sdss/rp.pb', format='ascii', data_start=1)
ifilt = Table.read(directory+'sdss/ip.pb', format='ascii', data_start=1)


# Get correct indices for star and filter spectra
test = Table.read('../PICKLES/a0i.sed', format='ascii')
nuv_starind, nuv_filtind = [], []
B_starind, B_filtind = [], []
V_starind, V_filtind = [], []
U_starind, U_filtind = [], []
R_starind, R_filtind = [], []
I_starind, I_filtind = [], []

J_starind, J_filtind = [], []
H_starind, H_filtind = [], []
K_starind, K_filtind = [], []

u_starind, u_filtind = [], []
g_starind, g_filtind = [], []
r_starind, r_filtind = [], []
i_starind, i_filtind = [], []
inda, indn = [], []

for ii in range(len(test)):
    for iii in range(len(nuvfilt)):
        if test['col1'][ii] == nuvfilt['col1'][iii]:
            nuv_starind.append(ii)
            nuv_filtind.append(iii)

    for iii in range(len(Bfilt)):
        if test['col1'][ii] == Bfilt['col1'][iii]:
            B_starind.append(ii)
            B_filtind.append(iii)
    for iii in range(len(Vfilt)):
        if test['col1'][ii] == Vfilt['col1'][iii]:
            V_starind.append(ii)
            V_filtind.append(iii)
    for iii in range(len(Ufilt)):
        if test['col1'][ii] == Ufilt['col1'][iii]:
            U_starind.append(ii)
            U_filtind.append(iii)
    for iii in range(len(Rfilt)):
        if test['col1'][ii] == Rfilt['col1'][iii]:
            R_starind.append(ii)
            R_filtind.append(iii)
    for iii in range(len(Ifilt)):
        if test['col1'][ii] == Ifilt['col1'][iii]:
            I_starind.append(ii)
            I_filtind.append(iii)

    for iii in range(len(Jfilt)):
        if test['col1'][ii] == Jfilt['col1'][iii]:
            J_starind.append(ii)
            J_filtind.append(iii)
    for iii in range(len(Hfilt)):
        if test['col1'][ii] == Hfilt['col1'][iii]:
            H_starind.append(ii)
            H_filtind.append(iii)
    for iii in range(len(Kfilt)):
        if test['col1'][ii] == Kfilt['col1'][iii]:
            K_starind.append(ii)
            K_filtind.append(iii)

    for iii in range(len(ufilt)):
        if test['col1'][ii] == ufilt['col1'][iii]:
            u_starind.append(ii)
            u_filtind.append(iii)
    for iii in range(len(gfilt)):
        if test['col1'][ii] == gfilt['col1'][iii]:
            g_starind.append(ii)
            g_filtind.append(iii)
    for iii in range(len(rfilt)):
        if test['col1'][ii] == rfilt['col1'][iii]:
            r_starind.append(ii)
            r_filtind.append(iii)
    for iii in range(len(ifilt)):
        if test['col1'][ii] == ifilt['col1'][iii]:
            i_starind.append(ii)
            i_filtind.append(iii)

# Extract correct spectra at the correct wavelengths for each band
nuvfilt = nuvfilt[nuv_filtind]
Bfilt = Bfilt[B_filtind]
Vfilt = Vfilt[V_filtind]
Ufilt = Ufilt[U_filtind]
Rfilt = Rfilt[R_filtind]
Ifilt = Ifilt[I_filtind]

Jfilt = Jfilt[J_filtind]
Hfilt = Hfilt[H_filtind]
Kfilt = Kfilt[K_filtind]

ufilt = ufilt[u_filtind]
gfilt = gfilt[g_filtind]
rfilt = rfilt[r_filtind]
ifilt = ifilt[i_filtind]

# Compute sums for denominator: wavelength * S(wavelength)_filter
nuvfiltsum = np.sum(nuvfilt['col2'] / nuvfilt['col1']) * 3631 * c_light
Bfiltsum = np.sum(Bfilt['col2'] / Bfilt['col1']) * 3631 * c_light
Vfiltsum = np.sum(Vfilt['col2'] / Vfilt['col1']) * 3631 * c_light
Ufiltsum = np.sum(Ufilt['col2'] / Ufilt['col1']) * 3631 * c_light
Rfiltsum = np.sum(Rfilt['col2'] / Rfilt['col1']) * 3631 * c_light
Ifiltsum = np.sum(Ifilt['col2'] / Ifilt['col1']) * 3631 * c_light

Jfiltsum = np.sum(Jfilt['col2'])  # / Jfilt['col1'])
Hfiltsum = np.sum(Hfilt['col2'])  # / Hfilt['col1'])
Kfiltsum = np.sum(Kfilt['col2'])  # / Kfilt['col1'])

ufiltsum = np.sum(ufilt['col2'] / ufilt['col1']) * 3631 * c_light
gfiltsum = np.sum(gfilt['col2'] / gfilt['col1']) * 3631 * c_light
rfiltsum = np.sum(rfilt['col2'] / rfilt['col1']) * 3631 * c_light
ifiltsum = np.sum(ifilt['col2'] / ifilt['col1']) * 3631 * c_light

name, nuv, B, V, U, R, IB, J, H, K, u, g, r, ib = [], [], [], [], [], [], [], [], [], [], [], [], [], []
picklefiles = np.loadtxt('../PICKLES/pickles.txt', dtype='string')

# Compute numerator: wavelength * S()
for i in picklefiles:
    table = Table.read('../PICKLES/'+i, format='ascii')
    name.append(i[:-4])
    nuv.append(np.sum(table['col2'][nuv_starind] * nuvfilt['col2'] * nuvfilt['col1']))
    B.append(np.sum(table['col2'][B_starind] * Bfilt['col2'] * Bfilt['col1']))
    V.append(np.sum(table['col2'][V_starind] * Vfilt['col2'] * Vfilt['col1']))
    U.append(np.sum(table['col2'][U_starind] * Ufilt['col2'] * Ufilt['col1']))
    R.append(np.sum(table['col2'][R_starind] * Rfilt['col2'] * Rfilt['col1']))
    IB.append(np.sum(table['col2'][I_starind] * Ifilt['col2'] * Ifilt['col1']))

    #j.append(np.sum(table['col2'][J_starind] * Jfilt['col2'] * Jfilt['col1']))
    #h.append(np.sum(table['col2'][H_starind] * Hfilt['col2'] * Hfilt['col1']))
    #k.append(np.sum(table['col2'][K_starind] * Kfilt['col2'] * Kfilt['col1']))
    J.append(np.sum(table['col2'][J_starind] * Jfilt['col2']))
    H.append(np.sum(table['col2'][H_starind] * Hfilt['col2']))
    K.append(np.sum(table['col2'][K_starind] * Kfilt['col2']))

    u.append(np.sum(table['col2'][u_starind] * ufilt['col2'] * ufilt['col1']))
    g.append(np.sum(table['col2'][g_starind] * gfilt['col2'] * gfilt['col1']))
    r.append(np.sum(table['col2'][r_starind] * rfilt['col2'] * rfilt['col1']))
    ib.append(np.sum(table['col2'][i_starind] * ifilt['col2'] * ifilt['col1']))


# Now compute f_lambda
nuv = [x/nuvfiltsum for x in nuv]
B = [x/Bfiltsum for x in B]
V = [x/Vfiltsum for x in V]
U = [x/Ufiltsum for x in U]
R = [x/Rfiltsum for x in R]
IB = [x/Ifiltsum for x in IB]

J = [x/Jfiltsum for x in J]
H = [x/Hfiltsum for x in H]
K = [x/Kfiltsum for x in K]

u = [x/ufiltsum for x in u]
g = [x/gfiltsum for x in g]
r = [x/rfiltsum for x in r]
ib = [x/ifiltsum for x in ib]

pbands = Table([nuv, U, B, V, R, IB, J, H, K, u, g, r, ib], names=('nuv', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'u', 'g', 'r', 'i'))

'''
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
Blamb = np.sqrt(np.sum(Bfilt['col2'])/(np.sum(Bfilt['col2']*(Bfilt['col1']+0.0)**-2)))
Vlamb = np.sqrt(np.sum(Vfilt['col2'])/(np.sum(Vfilt['col2']*(Vfilt['col1']+0.0)**-2)))
Ulamb = np.sqrt(np.sum(Ufilt['col2'])/(np.sum(Ufilt['col2']*(Ufilt['col1']+0.0)**-2)))
Rlamb = np.sqrt(np.sum(Rfilt['col2'])/(np.sum(Rfilt['col2']*(Rfilt['col1']+0.0)**-2)))
Ilamb = np.sqrt(np.sum(Ifilt['col2'])/(np.sum(Ifilt['col2']*(Ifilt['col1']+0.0)**-2)))

Jlamb = np.sqrt(np.sum(Jfilt['col2'])/(np.sum(Jfilt['col2']*(Jfilt['col1']+0.0)**-2)))
Hlamb = np.sqrt(np.sum(Hfilt['col2'])/(np.sum(Hfilt['col2']*(Hfilt['col1']+0.0)**-2)))
Klamb = np.sqrt(np.sum(Kfilt['col2'])/(np.sum(Kfilt['col2']*(Kfilt['col1']+0.0)**-2)))

ulamb = np.sqrt(np.sum(ufilt['col2'])/(np.sum(ufilt['col2']*(ufilt['col1']+0.0)**-2)))
glamb = np.sqrt(np.sum(gfilt['col2'])/(np.sum(gfilt['col2']*(gfilt['col1']+0.0)**-2)))
rlamb = np.sqrt(np.sum(rfilt['col2'])/(np.sum(rfilt['col2']*(rfilt['col1']+0.0)**-2)))
ilamb = np.sqrt(np.sum(ifilt['col2'])/(np.sum(ifilt['col2']*(ifilt['col1']+0.0)**-2)))


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
B = -2.5*np.log10(pbands['B'])
V = -2.5*np.log10(pbands['V'])
U = -2.5*np.log10(pbands['U'])
R = -2.5*np.log10(pbands['R'])
IB = -2.5*np.log10(pbands['I'])

J = -2.5*np.log10(pbands['J'] * (Jlamb)**2) + 7.96
H = -2.5*np.log10(pbands['H'] * (Hlamb)**2) + 7.526
K = -2.5*np.log10(pbands['K'] * (Klamb)**2) + 7.060
#j = -2.5*np.log10(pbands['j']) + 7.96
#h = -2.5*np.log10(pbands['h']) + 7.526
#k = -2.5*np.log10(pbands['k']) + 7.060

u = -2.5*np.log10(pbands['u'])
g = -2.5*np.log10(pbands['g'])
r = -2.5*np.log10(pbands['r'])
ib = -2.5*np.log10(pbands['i'])


data = [name, nuv, U, B, V, R, IB, J, H, K, u, g, r, ib]
ascii.write(data, '../picklemags_laphare_final.txt', names=['name', 'nuv', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'u', 'g', 'r', 'i'])
