from astropy.io import ascii
import numpy as np
from astropy.table import Table

pbands = Table.read('../picklesbands.txt', format='ascii')
name = pbands['name']

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
ascii.write(data, 'picklemags.txt', names=['name', 'nuv', 'b', 'v', 'j', 'h', 'k'])
