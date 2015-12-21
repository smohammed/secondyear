#WD distances
import numpy as np
from astropy.table import Table

v = Table.read('wds_vphasonly.txt',format='ascii')
wd = Table.read('wds_gais_vphas_newgrcut.txt',format='ascii')
vcut = np.where((v['gl'] > 200) & (v['gl'] < 250))
wdcut = np.where((wd['gl_galex'] > 200) & (wd['gl_galex'] < 250))

v = v[vcut]
wd = wd[wdcut]

sigma = 5.6704 * 10**-5     # erg/cm^2/s/K^4
radius = 637.1 * 10**6/2      # cm
Temp = 30 * 10**3           # K
Lsun = 3.846*10**33         # erg/s
Lum = sigma * Temp**4 * radius**2 * 4 * np.pi   # erg/s

Mabs = -2.5 * np.log10(Lum/Lsun)

DL1 = 10**((wd['nuv_mag']-Mabs)/5. + 1)  # parsecs

freq = 3*10**10/(2.267*10**-5)
DL2 = np.sqrt(Lum/(4*np.pi* 10**(wd['nuv_mag']/-2.5)*3631.*10**-23 * (freq))) / (3.086*10**18)  # parsecs

