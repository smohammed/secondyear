import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import hstack

#glinit, glend = 10, 20
#gbinit, gbend = 0, 10

###########################################################
# Takes a box, gets the tycho stars in it, takes photometry, returns residuals with gl,gb positions
###########################################################
def phot(glinit, glend, gbinit, gbend):
    intmap = fits.open('../combmaps12000/comb12000_gl'+str(glinit)+'to'+str(glend)+'.fits')[0].data
    tycho = fits.open('../tychobstarmatch.fits')[1].data
    tcut = np.where((tycho.gl > glinit) & (tycho.gl < glend) & (tycho.gb > gbinit) & (tycho.gb < gbend) & (tycho.nuv_cps > 5.) & (tycho.nuv_cps < 100.))
    tycho = tycho[tcut]
    tychogl = (tycho.gl - glinit) * 12000. / 20.
    tychogb = (tycho.gb - gbinit) * 12000. / 20.

    # Aperture radii
    r_search = 10
    r_in = 12
    r_out = 15

    #positions = [(tychogl,tychogb)]
    positions = []
    for i in range(len(tychogl)):
        positions.append((tychogl[i], tychogb[i]))

    apertures = CircularAperture(positions, r=r_search)

    # Background subtraction
    # Create annuli for subtraction range
    annulus_apertures = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

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
    #print(phot_table['residual_aperture_sum'])

    #plt.scatter(np.log10(phot_table['residual_aperture_sum']), np.log10(tycho.nuv_cps),edgecolor='none',alpha=0.5)
    #plt.xlabel('log10(Residual aperture sum)')
    #plt.ylabel('log10(nuv_cps)')
    #plt.title('gl =('+str(glinit)+','+str(glend)+'),gb =('+str(gbinit)+','+str(gbend)+'), r_search = 10 ,r_bkgd = 12 to 15')
    #plt.savefig('../gl =('+str(glinit)+','+str(glend)+'),gb =('+str(gbinit)+','+str(gbend)+').png')
    #plt.show()

    return tycho.gl, tycho.gb, phot_table['residual_aperture_sum']
