from astropy.io import fits, ascii
import numpy as np
from astropy.table import Table
import scipy.ndimage as ndimage
import scipy.stats as stats

# Total scans
#scans = ['0-10', '18-28', '27-37', '36-46', '45-55', '63-73', '108-118', '117-127', 

scans = ['126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '351-1', '315-325', '324-334', '333-343', '342-352', '72-82', '81-91', '90-100', '9-19', '99-109']


imeans = []
cmeans = []
smoothed = 0

for scan in scans:
    print('Running '+scan)
    # Get cmed and imed
    hdu = fits.open('../../galexscans/count_map_'+scan+'_in.fits')
    img = hdu[0].data

    if smoothed == 1:
        del img
        img = fits.open('../../galexscans/im1_'+scan+'.fits')[0].data

    ct = fits.open('../../galexscans/count_map_'+scan+'_count.fits')[0].data

    struc3 = ndimage.generate_binary_structure(2,2).astype(img.dtype)
    mask = ndimage.binary_dilation(img, structure=struc3, iterations=1)

    im1 = img[mask]
    ct1 = ct[mask]

    statimg = stats.sigmaclip(im1, 3, 3)
    statct = stats.sigmaclip(ct1, 3, 3)
    
    del img
    del ct
    del ct1
    del im1

    imean = np.mean(statimg[0])
    cmean = np.mean(statct[0])

    print('imean =', imean)
    print('cmean =', cmean)

   
    imeans.append(imean)
    cmeans.append(cmean)

maskvals = Table([scans, imeans, cmeans], names=['scan', 'imean', 'cmean'])

if smoothed == 0:
    ascii.write(maskvals, '../../galexscans/mean_scanvals.txt', format='basic')

if smoothed == 1:
    ascii.write(maskvals, '../../galexscans/mean_scanvals_smoothed.txt', format='basic')

print imeans
print cmeans
