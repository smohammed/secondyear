from astropy.io import fits
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy import stats
import os
matplotlib.rcParams['figure.figsize'] = 8,8

path = "../corrcsv/"
files = np.loadtxt("../filelists/csv_galcorr.txt", dtype = 'string')

for line in range(len(files)):
	galscan = fits.open(path+files[line])[1].data

	time = galscan.time
	uniquetime = np.unique(np.round(np.unique(time)))

	for i in uniquetime:
		gallim = np.where((time >= uniquetime[i]) & (time < uniquetime[i+1]))

		gl = galscan.gl[gallim]
		gb = galscan.gb[gallim]

		for step in range(len(gl)):
			if min(gl) < 10. and gl[step] > 310.:
				gl[step] = gl[step] - 360.

		xmin,xmax = min(gl),max(gl)
		ymin,ymax = min(gb),max(gb)
		X,Y = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
		
		positions = np.vstack([X.ravel(),Y.ravel()])
		values = np.vstack([gl,gb])
		kernel = stats.gaussian_kde(values)
		Z = np.reshape(kernel(positions).T,X.shape)

		#plt.plot(gl,gb,'k.',markersize=1)
		plt.imshow(np.rot90(Z),extent=[-0.5,1.5,ymin,ymax])
		
		if not os.path.exists('../kdemaps/'+str(files[line][:5])+'corrcsv/'):
			os.makedirs('../kdemaps/'+str(files[line][:5])+'corrcsv/')

		#plt.show()
		plt.savefig('../kdemaps/'+str(files[line][:5])+'corrcsv/'+str(files[line][:5])+'gaus_kde_'+str(i)+'.png')
		print str(files[line][:5])+' for time '+str(i)



'''
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

path = "../corrcsv/"
files = np.loadtxt("../filelists/csv_galcorr.txt", dtype = 'string')
a = fits.open(path+files[0])[1].data
alim = np.where((a.time >= 30.) & (a.time < 31.))

gl = a.gl[alim]
gb = a.gb[alim]
time = a.time[alim]
for i in range(len(gl)):
	if min(gl) < 10. and gl[i] > 310.:
		gl[i] = gl[i] - 360.

xmin,xmax = min(gl),max(gl)
ymin,ymax = min(gb),max(gb)
X,Y = np.mgrid[xmin:xmax:100j,ymin:ymax:100j]
positions = np.vstack([X.ravel(),Y.ravel()])
values = np.vstack([gl,gb])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T,X.shape)
plt.imshow(np.rot90(Z),extent=[xmin,xmax,ymin,ymax])
plt.show()
'''