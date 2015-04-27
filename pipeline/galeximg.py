from astropy.io import fits
import numpy as np
from numpy import loadtxt
from astropy import wcs

########################################################################################
# This code generates an intensity map. Combine with the exposure map to get a final image.
########################################################################################

############################################
#Generate which pointing map?
############################################
pointing = 2
correction = 0

############################################
# Input galactic longitude range
############################################
a = 10				# Start gl
b = 20				# End gl
c = 0	 			# Start gb
d = 10 				# End gb

def intensitymap(glinit,glfinal,gbinit,gbfinal,pointing,correction):
	
	############################################
	# Set the file pathways
	############################################
	
	if pointing == 1:
		path = "../intmapfiles/"
		lines = loadtxt("../filelists/pointingtimes.txt", comments = "#", delimiter = ",", unpack = False,	 dtype = str)
	
	if pointing == 2:
		path = "../corrcsv/"
		#lines = loadtxt("../filelists/AIS_GAL_SCAN_csv.txt", comments = "#", delimiter = ",", unpack = 	False, dtype = str)
		lines = loadtxt("../filelists/csv_galcorr.txt", dtype = 'string')
	
	if correction == 1:
		offpath = "../dxdyoffsets/"
		offlines = loadtxt("../filelists/dxdycorrectionfiles.txt",dtype='string')
	
	############################################
	# Grab appropriate files using gl range
	############################################
	startgl = 0
	endgl = 0
	
	if pointing == 1:
		while (float(lines[startgl].split('_')[0])/10.) < a:
			startgl+=1
		while (float(lines[endgl].split('_')[0])/10.) < b:
			endgl+=1
		endgl = endgl - 1
		
		photons = fits.open(path+lines[startgl])[1].data
		photons2 = fits.open(path+lines[startgl+1])[1].data
		
		gl = np.concatenate((photons.gl_cor,photons2.gl_cor))
		gb = np.concatenate((photons.gb_cor,photons2.gb_cor))
	
	if pointing == 2:
		while (float(lines[startgl].split('_')[0])/10.) < a:
			startgl+=1
		while (float(lines[endgl].split('_')[0])/10.) < b:
			endgl+=1
		endgl = endgl - 1
		
		photons = fits.open(path+lines[startgl])[1].data
		photlim = np.where((photons.gl >= a) & (photons.gl < b) & (photons.gb >= c) & (photons.gb < d))
		gl = photons.gl[photlim]
		gb = photons.gb[photlim]
		time = photons.time[photlim]
	
		gl[(min(gl) < 10.) & (gl > 350.)] = gl[(min(gl) < 10.) & (gl > 350.)] - 360.
	
		if correction == 1:
			offset = fits.open(offpath+offlines[startgl])[1].data
	
			for i in range(len(offset)-1):
				offsetind = np.where((time >= offset.time[i]) & (time < offset.time[i+1]))
				gl[offsetind] = gl[offsetind] - offset.dx[i]/60.
				gb[offsetind] = gb[offsetind] - offset.dy[i]/60.
	
	
	print 'startgl = ', startgl
	print 'endgl =', endgl
	
	############################################
	# Combine photons from up and down runs, stack data horizontally
	############################################
	if pointing == 1:
		for i in range(startgl+2,endgl,2):
			print i
			try:
				f1 = fits.open(path+lines[i])[1].data
				f1gl_cor = f1.gl_cor
				f1gb_cor = f1.gb_cor
		
			except IOError:								# if the file is corrupted, skip it 
				print 'Index '+str(i)+' is corrupted'
		
			try:
				f2 = fits.open(path+lines[i+1])[1].data
				f2gl_cor = f2.gl_cor
				f2gb_cor = f2.gb_cor
		
			except IOError:								# if the file is corrupted, skip it 
				print 'Index '+str(i+1)+' is corrupted'
			
			glinit = np.concatenate((f1gl_cor,f2gl_cor))
			gbinit = np.concatenate((f1gb_cor,f2gb_cor))
			gl = np.hstack([gl, glinit])
			gb = np.hstack([gb, gbinit])
	
	if pointing == 2:
		for i in range(startgl+1,endgl+1):
			print i
			try:
				f1 = fits.open(path+lines[i])[1].data 
				f1lim = np.where((f1.gl >= a) & (f1.gl < b) & (f1.gb >= c) &  (f1.gb < d))
				f1gl = f1.gl[f1lim]
				f1gb = f1.gb[f1lim]
				time = f1.time[f1lim]
				f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] = f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] - 360.
	
				if correction == 1:
					offset= fits.open(offpath+offlines[i])[1].data
					for i in range(len(offset)-1):
						offsetind = np.where((time >= offset.time[i]) & (time < offset.time[i+1]))
						f1gl[offsetind] = f1gl[offsetind] - offset.dx[i]/60.
						f1gb[offsetind] = f1gb[offsetind] - offset.dy[i]/60.
	
			except IOError:								# if the file is corrupted, skip it 
				print 'Index '+str(i)+' is corrupted'
				
			gl = np.hstack([gl, f1gl])
			gb = np.hstack([gb, f1gb])
	
	############################################
	# Bin data and plot stuff
	############################################
	if pointing == 1:
		binnum = 1200
		
		H, xbins, ybins = np.histogram2d(gl, gb, bins = (np.linspace(a, b, binnum), np.linspace(c, d, 	binnum)))
		#fig = plt.figure(figsize = (10,10))
		#ax = plt.axes()
		
		fits.PrimaryHDU(H).writeto('../intensitymap'+str(binnum)+'_'+str(a)+'_'+str(b)+'.fits')
		
		'''
		ax.imshow(np.sqrt(H).T, vmin = 0, vmax = 5, origin = 'lower', extent = [a, b, -10, 10], 	interpolation = 'nearest', aspect = 'auto', cmap = 'gray')
		ax.set_xlabel(r'${\rm gl}$')
		ax.set_ylabel(r'${\rm gb}$')
		ax.set_title('Galactic Plane in UV, %d<gl<%d'%(a,b))
		
		fig.set_dpi(300)
		fig.savefig('GALEX_galplane_%d-%d_vmax05.png'%(a,b), bbox_inches = 'tight', pad_inches = 0)
		''' 
	
	if pointing == 2:
		binnum = 12000
		print binnum 
	
		H, xbins, ybins = np.histogram2d(gl, gb, bins = (np.linspace(a, b, binnum), np.linspace(c, d, 	binnum)))
		print np.size(H)
	
	
	
		if correction == 0:
			return fits.PrimaryHDU(H.T).writeto('../intmapcsvcorr'+str(binnum)+'_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits')
		elif correction == 1:
			return fits.PrimaryHDU(H).writeto('../intmapcsvcorr'+str(binnum)+'_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'_corr_no_offset.fits')
		
		'''
		fig = plt.figure(figsize = (10,10))
		ax = plt.axes()
		ax.imshow(np.sqrt(H).T, vmin = 0, vmax = 5, origin = 'lower', extent = [a, b, -10, 10], 	interpolation = 'nearest', aspect = 'auto', cmap = 'gray')
		ax.set_xlabel(r'${\rm gl}$')
		ax.set_ylabel(r'${\rm gb}$')
		ax.set_title('Galactic Plane in UV, %d<gl<%d'%(a,b))
		
		fig.set_dpi(300)
		fig.savefig('GALEX_galplane_%d-%d_vmax05.png'%(a,b), bbox_inches = 'tight', pad_inches = 0)
		'''
	
for glstep in range(30, 361, 10):
	intensitymap(glstep,glstep+10,0,10,pointing,correction)