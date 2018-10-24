from astropy.io import fits
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

scans = ['9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']

#################################################
# Load catalogs
#################################################
gg = fits.open('../plane_gais_09-26-18.fits')[1].data
dx = gg['glon']-gg['gl']
dy = gg['glat']-gg['gb']

#dx[np.where(dx < -100.)] = dx[np.where(dx < -100.)] + 360

# Label partial or full scans
for region in range(0, 360, 5):
    print region
    q = np.where((gg['glon'] > region) & (gg['glon'] < (region + 5)))
    dxq = dx[q]
    dyq = dy[q]

    plt.scatter(dxq*3600., dyq*3600., s=1, marker='.')

    plt.axvline(x=0, color='black')
    plt.axhline(y=0,  color='black')

    plt.xlabel('dx [GAIS - Plane] [arcsec]')
    plt.ylabel('dy [GAIS - Plane] [arcsec]')
    plt.xlim(-3.5, 3.5)
    plt.ylim(-3.5, 3.5)
    plt.title(str(region)+'-'+str(region+5))

    plt.savefig('../images/gaisplaneoffsets/09-26-18_offsets_'+str(region)+'-'+str(region+5)+'.png')

    #plt.show()
    plt.close()
