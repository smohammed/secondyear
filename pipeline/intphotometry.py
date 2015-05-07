import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
#from astropy import units as u
#from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16,8
matplotlib.rcParams['font.size'] = 20
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import hstack

glinit,glend = 10,20
gbinit,gbend = 0,10
intmap = np.sqrt(fits.open('../intmapcsvcorr12000_gl_'+str(glinit)+'to'+str(glend)+'_gb_'+str(gbinit)+'to'+str(gbend)+'transp.fits')[0].data)
t2 = fits.open('../tychobstarmatch.fits')[1].data
tcut = np.where((t2.gl > glinit) & (t2.gl < glend) & (t2.gb > gbinit) & (t2.gb < gbend))
t2 = t2[tcut]
t2gl = (t2.gl - glinit) * 12000. / 10.
t2gb = (t2.gb - gbinit) * 12000. / 10.


# Photometry
r_search = 10
r_in = 12
r_out = 15

#positions = [(t2gl,t2gb)]
positions = []
for i in range(len(t2gl)):
	positions.append((t2gl[i],t2gb[i]))

apertures = CircularAperture(positions, r=10.)

# Background subtraction
# Create annuli for subtraction range
annulus_apertures = CircularAnnulus(positions, r_in=12., r_out=15.)

# Compute the sum inside the inner aperture and then the bkgd flux
rawflux_table = aperture_photometry(intmap, apertures)
bkgflux_table = aperture_photometry(intmap, annulus_apertures)
phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])

# calculate the mean local background
bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()

# Finally subtract background
bkg_sum = bkg_mean * apertures.area()
final_sum = phot_table['aperture_sum_raw'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum
print(phot_table['residual_aperture_sum'])

#plt.scatter(np.log10(phot_table['residual_aperture_sum']), np.log10(t2.nuv_cps),edgecolor='none',alpha=0.5)
#plt.xlabel('log10(Residual aperture sum)')
#plt.ylabel('log10(nuv_cps)')
#plt.title('gl =('+str(glinit)+','+str(glend)+'),gb =('+str(gbinit)+','+str(gbend)+'), r_search = 10 ,r_bkgd = 12 to 15')
#plt.savefig('../gl =('+str(glinit)+','+str(glend)+'),gb =('+str(gbinit)+','+str(gbend)+').png')
#plt.show()

