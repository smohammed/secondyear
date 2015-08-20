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

# Open LePHARE files
directory = '../PICKLES/filt/'
nuvfilt = Table.read(directory+'galex/NUV.pb', format='ascii')
bfilt = Table.read(directory+'wfi/B.pb', format='ascii', data_start=1)
vfilt = Table.read(directory+'wfi/V.pb', format='ascii', data_start=1)
jfilt = Table.read(directory+'2mass/J.pb', format='ascii', data_start=1)
hfilt = Table.read(directory+'2mass/H.pb', format='ascii', data_start=1)
kfilt = Table.read(directory+'2mass/Ks.pb', format='ascii', data_start=1)


# Get correct indices
test = Table.read('../PICKLES/a0i.sed', format='ascii')
nuv_starind, nuv_filtind = [], []
b_starind, b_filtind = [], []
v_starind, v_filtind = [], []
j_starind, j_filtind = [], []
h_starind, h_filtind = [], []
k_starind, k_filtind = [], []
inda, indn = [], []

for i in range(len(test)):
    for j in range(len(nuvfilt)):
        if test['col1'][i] == nuvfilt['col1'][j]:
            nuv_starind.append(i)
            nuv_filtind.append(j)

    for j in range(len(bfilt)):
        if test['col1'][i] == bfilt['col1'][j]:
            b_starind.append(i)
            b_filtind.append(j)

    for j in range(len(vfilt)):
        if test['col1'][i] == vfilt['col1'][j]:
            v_starind.append(i)
            v_filtind.append(j)

    for j in range(len(jfilt)):
        if test['col1'][i] == jfilt['col1'][j]:
            j_starind.append(i)
            j_filtind.append(j)

    for j in range(len(hfilt)):
        if test['col1'][i] == hfilt['col1'][j]:
            h_starind.append(i)
            h_filtind.append(j)

    for j in range(len(kfilt)):
        if test['col1'][i] == kfilt['col1'][j]:
            k_starind.append(i)
            k_filtind.append(j)

nuvfilt = nuvfilt[nuv_filtind]
bfilt = bfilt[b_filtind]
vfilt = vfilt[v_filtind]
jfilt = jfilt[j_filtind]
hfilt = hfilt[h_filtind]
kfilt = kfilt[k_filtind]

# Compute sums for denominator
nuvfiltsum = np.sum(nuvfilt['col2'])
bfiltsum = np.sum(bfilt['col2'])
vfiltsum = np.sum(vfilt['col2'])
jfiltsum = np.sum(jfilt['col2'])
hfiltsum = np.sum(hfilt['col2'])
kfiltsum = np.sum(kfilt['col2'])


# Now compute numerator
name, nuv, b, v, j, h, k = [], [], [], [], [], [], []
picklefiles = np.loadtxt('../PICKLES/pickles.txt', dtype='string')

for i in picklefiles:
    table = Table.read('../PICKLES/'+i, format='ascii')
    name.append(i[:-4])
    nuv.append(np.sum(table['col2'][nuv_starind] * nuvfilt['col2']))
    b.append(np.sum(table['col2'][b_starind] * bfilt['col2']))
    v.append(np.sum(table['col2'][v_starind] * vfilt['col2']))
    j.append(np.sum(table['col2'][j_starind] * jfilt['col2']))
    h.append(np.sum(table['col2'][h_starind] * hfilt['col2']))
    k.append(np.sum(table['col2'][k_starind] * kfilt['col2']))

# Now compute f_lambda
nuv = [x/nuvfiltsum for x in nuv]
b = [x/bfiltsum for x in b]
v = [x/vfiltsum for x in v]
j = [x/jfiltsum for x in j]
h = [x/hfiltsum for x in h]
k = [x/kfiltsum for x in k]

# Old way
'''
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
'''

pbands = Table([nuv, b, v, j, k, k], names=('nuv', 'b', 'v', 'j', 'h', 'k'))

##########################################################################
# convert f_lambda to f_v, divide by 0 point, then conv to Jy. 1 Jy = 10^-23 erg/s/Hz/cm^2
##########################################################################

#nuv = -2.5*np.log10((((pbands['nuv']-(2.22*10**-8))*3.34*10**4*(2300)**2))*10**-23)-48.6

# Wavelengths in A
#nuv = -2.5*np.log10(((pbands['nuv'] * 3.34*10**4 * (2300)**2)) * 10**-23) - 48.6
#b = -2.5*np.log10(pbands['b'] * 3.34*10**4 * (4350)**2 * 10**-23) - 48.6
#v = -2.5*np.log10(pbands['v'] * 3.34*10**4 * (5050)**2 * 10**-23) - 48.6
#j = -2.5*np.log10(((pbands['j']*3.34*10**4 * (12350)**2) / 1594) * 10**-23) - 48.6
#h = -2.5*np.log10(((pbands['h']*3.34*10**4 * (16620)**2) / 1024) * 10**-23) - 48.6
#k = -2.5*np.log10(((pbands['k']*3.34*10**4 * (21590)**2) / 666.7) * 10**-23) - 48.6

nuv = -2.5*np.log10(((pbands['nuv'] * 3.34*10**4 * (2300)**2))) + 8.9

b = -2.5*np.log10(pbands['b'] * 3.34*10**4 * (4350)**2) + 8.9
v = -2.5*np.log10(pbands['v'] * 3.34*10**4 * (5050)**2) + 8.9

j = -2.5*np.log10(((pbands['j']*3.34*10**4 * (12350)**2))) + 7.96
h = -2.5*np.log10(((pbands['h']*3.34*10**4 * (16620)**2))) + 7.526
k = -2.5*np.log10(((pbands['k']*3.34*10**4 * (21590)**2))) + 7.060


data = [name, nuv, b, v, j, h, k]
ascii.write(data, '../picklemags_laphare.txt', names=['name', 'nuv', 'b', 'v', 'j', 'h', 'k'])
