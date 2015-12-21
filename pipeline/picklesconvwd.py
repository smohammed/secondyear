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
# col 1 = angstroms, col2 = erg/s/cm^2/A
directory = '../PICKLES/filt/'
nuvfilt = Table.read(directory+'galex/NUV.pb', format='ascii')
bfilt = Table.read(directory+'wfi/B.pb', format='ascii', data_start=1)
vfilt = Table.read(directory+'wfi/V.pb', format='ascii', data_start=1)

ufilt = Table.read(directory+'sdss/up.pb', format='ascii', data_start=1)
gfilt = Table.read(directory+'sdss/gp.pb', format='ascii', data_start=1)
rfilt = Table.read(directory+'sdss/rp.pb', format='ascii', data_start=1)
ifilt = Table.read(directory+'sdss/ip.pb', format='ascii', data_start=1)

name, nuv, b, v, u, g, r, ib = [], [], [], [], [], [], [], []

picklesfiles = np.loadtxt('../PICKLES/wds/pickleswdlist.txt', dtype='string')

for filename in picklesfiles:
    # Get correct indices
    test = Table.read('../PICKLES/wds/'+filename, format='fits')
    nuv_starind, nuv_filtind = [], []
    b_starind, b_filtind = [], []
    v_starind, v_filtind = [], []
    u_starind, u_filtind = [], []
    g_starind, g_filtind = [], []
    r_starind, r_filtind = [], []
    i_starind, i_filtind = [], []

    # Match wavelengths between filters and reference stars
    for i in range(len(test)):
        for j in range(len(nuvfilt)-1):
            if (test['WAVELENGTH'][i] > nuvfilt['col1'][j]) and (test['WAVELENGTH'][i] < nuvfilt['col1'][j+1]):
                nuv_starind.append(i)
                nuv_filtind.append(j)

        for j in range(len(bfilt)-1):
            if (test['WAVELENGTH'][i] > bfilt['col1'][j]) and (test['WAVELENGTH'][i] < bfilt['col1'][j+1]):
                b_starind.append(i)
                b_filtind.append(j)

        for j in range(len(vfilt)-1):
            if (test['WAVELENGTH'][i] > vfilt['col1'][j]) and (test['WAVELENGTH'][i] < vfilt['col1'][j+1]):
                v_starind.append(i)
                v_filtind.append(j)

        for j in range(len(ufilt)-1):
            if (test['WAVELENGTH'][i] > ufilt['col1'][j]) and (test['WAVELENGTH'][i] < ufilt['col1'][j+1]):
                u_starind.append(i)
                u_filtind.append(j)

        for j in range(len(gfilt)-1):
            if (test['WAVELENGTH'][i] > gfilt['col1'][j]) and (test['WAVELENGTH'][i] < gfilt['col1'][j+1]):
                g_starind.append(i)
                g_filtind.append(j)

        for j in range(len(rfilt)-1):
            if (test['WAVELENGTH'][i] > rfilt['col1'][j]) and (test['WAVELENGTH'][i] < rfilt['col1'][j+1]):
                r_starind.append(i)
                r_filtind.append(j)

        for j in range(len(ifilt)-1):
            if (test['WAVELENGTH'][i] > ifilt['col1'][j]) and (test['WAVELENGTH'][i] < ifilt['col1'][j+1]):
                i_starind.append(i)
                i_filtind.append(j)

    # Pull out indices
    nuvfilt = nuvfilt[nuv_filtind]
    bfilt = bfilt[b_filtind]
    vfilt = vfilt[v_filtind]
    ufilt = ufilt[u_filtind]
    gfilt = gfilt[g_filtind]
    rfilt = rfilt[r_filtind]
    ifilt = ifilt[i_filtind]

    # Make sure they're the same length
    print len(nuv_starind) - len(nuv_filtind)
    print len(b_starind) - len(b_filtind)
    print len(v_starind) - len(v_filtind)
    print len(u_starind) - len(u_filtind)
    print len(g_starind) - len(g_filtind)
    print len(r_starind) - len(r_filtind)
    print len(i_starind) - len(i_filtind)

    # Compute sums for denominator. Units = erg/s/cm^2/A^2
    nuvfiltsum = np.sum(nuvfilt['col2']/nuvfilt['col1'])
    bfiltsum = np.sum(bfilt['col2']/bfilt['col1'])
    vfiltsum = np.sum(vfilt['col2']/vfilt['col1'])
    ufiltsum = np.sum(ufilt['col2']/ufilt['col1'])
    gfiltsum = np.sum(gfilt['col2']/gfilt['col1'])
    rfiltsum = np.sum(rfilt['col2']/rfilt['col1'])
    ifiltsum = np.sum(ifilt['col2']/ifilt['col1'])

    # Now compute numerator
    name1, nuv1, b1, v1, u1, g1, r1, ib1 = [], [], [], [], [], [], [], []

    # Units = erg/s/cm^2/A * erg/s/cm^2/A * A
    name1.append(filename[:-4])
    nuv1.append(np.sum(test['FLUX'][nuv_starind] * nuvfilt['col2']*nuvfilt['col1']))
    b1.append(np.sum(test['FLUX'][b_starind] * bfilt['col2']*bfilt['col1']))
    v1.append(np.sum(test['FLUX'][v_starind] * vfilt['col2']*vfilt['col1']))
    u1.append(np.sum(test['FLUX'][u_starind] * ufilt['col2']*ufilt['col1']))
    g1.append(np.sum(test['FLUX'][g_starind] * gfilt['col2']*gfilt['col1']))
    r1.append(np.sum(test['FLUX'][r_starind] * rfilt['col2']*rfilt['col1']))
    ib1.append(np.sum(test['FLUX'][i_starind] * ifilt['col2']*ifilt['col1']))

    # Now compute f_lambda. Units = erg/s/cm^2/A^2
    nuv1 = [x/nuvfiltsum for x in nuv1]
    b1 = [x/bfiltsum for x in b1]
    v1 = [x/vfiltsum for x in v1]
    u1 = [x/ufiltsum for x in u1]
    g1 = [x/gfiltsum for x in g1]
    r1 = [x/rfiltsum for x in r1]
    ib1 = [x/ifiltsum for x in ib1]

    name.append(name1)
    nuv.append(nuv1)
    b.append(b1)
    v.append(v1)
    u.append(u1)
    g.append(g1)
    r.append(r1)
    ib.append(ib1)


print 'name =', len(name)
print 'nuv =', len(nuv)
print 'b =', len(b)
print 'v =', len(v)
print 'u =', len(u)
print 'g =', len(g)
print 'r =', len(r)
print 'i =', len(ib)

ib[np.isnan(ib)] = 0

# Now to combine all the data
pbands = Table([name, np.array(nuv), np.array(b), np.array(v), np.array(u), np.array(g), np.array(r), np.array(ib)], names=('name','nuv', 'b', 'v', 'u', 'g', 'r', 'i'))

##########################################################################
# convert f_lambda to f_v, divide by 0 point, then conv to Jy. 1 Jy = 10^-23 erg/s/Hz/cm^2
##########################################################################
pbands['nuv'] = -2.5*np.log10(pbands['nuv']) + 8.9
pbands['b'] = -2.5*np.log10(pbands['b']) + 8.9
pbands['v'] = -2.5*np.log10(pbands['v']) + 8.9
pbands['u'] = -2.5*np.log10(pbands['u']) + 8.9
pbands['g'] = -2.5*np.log10(pbands['g']) + 8.9
pbands['r'] = -2.5*np.log10(pbands['r']) + 8.9
pbands['i'] = -2.5*np.log10(pbands['i']) + 8.9


print pbands['nuv']
print pbands['b']
print pbands['v']
print pbands['u']
print pbands['g']
print pbands['r']
print pbands['i']

ascii.write(pbands, '../picklemags_wds.txt', format='basic')

a = -2.5*np.log10(ib) + 8.9

a = [a[0][0],a[1][0],a[2][0],a[3][0],a[4][0],a[5][0]]

pbands.remove_column('i')
pbands['i'] = a