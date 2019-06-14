from astropy.io import fits
import matplotlib
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from astropy.wcs import wcs
matplotlib.rcParams['figure.figsize'] = 14, 12
matplotlib.rcParams['font.size'] = 18

# Grab a snapshot of a given coordinate in the UVGAPS images


scans = ['0-10', '9-19', '18-28', '27-37', '36-46', '45-55', '63-73', '72-82', '81-91', '90-100', '99-109', '108-118', '117-127', '126-136', '135-145', '144-154', '153-163', '162-172', '171-181', '180-190', '189-199', '198-208', '207-217', '216-226', '225-235', '234-244', '243-253', '252-262', '261-271', '270-280', '279-289', '288-298', '297-307', '306-316', '315-325', '324-334', '333-343', '342-352', '351-1']


scans1 = [0, 9, 18, 27, 36, 45, 63, 72, 81, 90, 99 , 108, 117, 126, 135, 144, 153, 162, 171, 180, 189, 198, 207, 216, 225, 234, 243, 252, 261, 270, 279, 288, 297, 306, 315, 324, 333, 342, 351] 

scans2 = [10, 19, 28, 37, 46, 55, 73, 82, 91, 100, 109, 118, 127, 136, 145, 154, 163, 172, 181, 190, 199, 208, 217, 226, 235, 244, 253, 262, 271, 280, 289, 298, 307, 316, 325, 334, 343, 352, 1]

def snap(gl, gb):
	# Find out which is the best image to grab the snapshot from
	'''
	for reg in range(len(scans)-1):
		if scans[reg] > scans[reg+1] and scans[reg] != 351:
			continue
		#print reg
	'''
	for reg in range(len(scans1)):
		if gl > scans1[reg] and gl < scans2[reg]:
			scan = scans[reg]

	# Load the image
	hdulist = fits.open('../../galexscans/count_map_'+scan+'_in.fits')
	img = hdulist[0].data
	# Convert gl, gb to pixels and extract the image

	w = wcs.WCS(hdulist[0].header)
	world = np.array([[gl, gb]])
	pix = w.wcs_world2pix(world, 1)

	im = img[int(pix[0][1])-100:int(pix[0][1])+100, int(pix[0][0])-100:int(pix[0][0])+100]
	wcsmap = w[int(pix[0][1])-100:int(pix[0][1])+100, int(pix[0][0])-100:int(pix[0][0])+100]
	header = wcsmap.to_header()

	fig = plt.figure()
	fig.add_subplot(111, projection=wcsmap)
	plt.imshow(im, cmap=cm.gray, origin='lower', vmin=0, vmax=0.7)
	plt.xlabel('nuv = '+str(ob2['nuv'][i])[:4])
	plt.title('sourceid = '+str(ob2['source_id'][i]))
	#plt.show()
	plt.savefig('../images/blueobj/planeobj_'+str(gl)[:4]+'_'+str(gb)[:4]+'.png')
	plt.close()




# Add:
# label all stars around object

