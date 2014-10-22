from astropy.io import fits
import numpy as np
from numpy import loadtxt
#from matplotlib import pyplot as plt

########################################################################################
# This code finds the photon times used in int and creates a list to compare with the expmap data
########################################################################################

############################################
# Input galactic longitude range
############################################
a=0 		# Start gl
b=20		# End gl

############################################
# Set the file pathways
############################################
path = "../intmapfiles/"
lines = loadtxt("../filelists/fits_photon_list.txt", comments = "#", delimiter = ",", unpack = False, dtype = str)

############################################
# Grab appropriate files using gl range
############################################
startgl = a
endgl = a

while (float(lines[startgl].split('_')[0])/10.) < a:
	startgl+=1
while (float(lines[endgl].split('_')[0])/10.) < b:
	endgl+=1
endgl = endgl - 1

photons = fits.open(path+lines[startgl])[1].data
photons2 = fits.open(path+lines[startgl+1])[1].data

time = np.concatenate((photons.time,photons2.time))

time = np.unique(time.astype(int))

np.savetxt('times/'+lines[startgl].split('_')[0]+'photontimes.txt',time,newline=',')


############################################
# Combine photons from up and down runs, stack data horizontally
############################################
for i in range(startgl+2,endgl,2):
	print i
	try:
		f1 = fits.open(path+lines[i])[1].data
		f1time = f1.time

	except IOError:								# if the file is corrupted, skip it 
		print 'Index '+str(i)+' is corrupted'

	try:
		f2 = fits.open(path+lines[i+1])[1].data
		f2time = f2.time

	except IOError:								# if the file is corrupted, skip it 
		print 'Index '+str(i+1)+' is corrupted'

	timeinit = np.unique(np.concatenate((f1time,f2time)).astype(int))

	np.savetxt('times/'+lines[i].split('_')[0]+'photontimes.txt',time,newline=',')