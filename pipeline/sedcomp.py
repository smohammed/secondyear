from astropy.io import ascii
import numpy as np
from astropy.table import Table
import pysynphot

# Filter files
directory = '/PICKLES/filt/'
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

# Pysynphot Pickles files
o5v_file = os.path.join(os.environ['PYSYN_CDBS'], 'grid/pickles/dat_uvi/', 'pickles_1.fits')
o5v = pysynphot.FileSpectrum(o5v_file)

nuvcut = np.where((o5v.wave > nuvfilt['col1'][0]) & (o5v.wave < nuvfilt['col1'][-1]))
Bcut = np.where((o5v.wave > Bfilt['col1'][0]) & (o5v.wave < Bfilt['col1'][-1]))
Vcut = np.where((o5v.wave > Vfilt['col1'][0]) & (o5v.wave < Vfilt['col1'][-1]))
Ucut = np.where((o5v.wave > Ufilt['col1'][0]) & (o5v.wave < Ufilt['col1'][-1]))
Rcut = np.where((o5v.wave > Rfilt['col1'][0]) & (o5v.wave < Rfilt['col1'][-1]))
Icut = np.where((o5v.wave > Ifilt['col1'][0]) & (o5v.wave < Ifilt['col1'][-1]))

Jcut = np.where((o5v.wave > Jfilt['col1'][0]) & (o5v.wave < Jfilt['col1'][-1]))
Hcut = np.where((o5v.wave > Hfilt['col1'][0]) & (o5v.wave < Hfilt['col1'][-1]))
Kcut = np.where((o5v.wave > Kfilt['col1'][0]) & (o5v.wave < Kfilt['col1'][-1]))

ucut = np.where((o5v.wave > ufilt['col1'][0]) & (o5v.wave < ufilt['col1'][-1]))
gcut = np.where((o5v.wave > gfilt['col1'][0]) & (o5v.wave < gfilt['col1'][-1]))
rcut = np.where((o5v.wave > rfilt['col1'][0]) & (o5v.wave < rfilt['col1'][-1]))
icut = np.where((o5v.wave > ifilt['col1'][0]) & (o5v.wave < ifilt['col1'][-1]))

nuvmag = -2.5*np.log10(np.sum(o5v.flux[nuvcut]*o5v.wave[nuvcut]*nuvfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[nuvcut]*nuvfilt['col2'])))
Bmag = -2.5*np.log10(np.sum(o5v.flux[Bcut]*o5v.wave[Bcut]*Bfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Bcut]*Bfilt['col2'])))
Vmag = -2.5*np.log10(np.sum(o5v.flux[Vcut]*o5v.wave[Vcut]*Vfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Vcut]*Vfilt['col2'])))
Umag = -2.5*np.log10(np.sum(o5v.flux[Ucut]*o5v.wave[Ucut]*Ufilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Ucut]*Ufilt['col2'])))
Rmag = -2.5*np.log10(np.sum(o5v.flux[Rcut]*o5v.wave[Rcut]*Rfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Rcut]*Rfilt['col2'])))
Imag = -2.5*np.log10(np.sum(o5v.flux[Icut]*o5v.wave[Icut]*Ifilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Icut]*Ifilt['col2'])))

Jmag = -2.5*np.log10(np.sum(o5v.flux[Jcut]*o5v.wave[Jcut]*Jfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Jcut]*Jfilt['col2'])))
Hmag = -2.5*np.log10(np.sum(o5v.flux[Hcut]*o5v.wave[Hcut]*Hfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Hcut]*Hfilt['col2'])))
Kmag = -2.5*np.log10(np.sum(o5v.flux[Kcut]*o5v.wave[Kcut]*Kfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[Kcut]*Kfilt['col2'])))

umag = -2.5*np.log10(np.sum(o5v.flux[ucut]*o5v.wave[ucut]*ufilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[ucut]*ufilt['col2'])))
gmag = -2.5*np.log10(np.sum(o5v.flux[gcut]*o5v.wave[gcut]*gfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[gcut]*gfilt['col2'])))
rmag = -2.5*np.log10(np.sum(o5v.flux[rcut]*o5v.wave[rcut]*rfilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[rcut]*rfilt['col2'])))
imag = -2.5*np.log10(np.sum(o5v.flux[icut]*o5v.wave[icut]*ifilt['col2'])/(np.sum(3631*3*10**10/o5v.wave[icut]*ifilt['col2'])))






