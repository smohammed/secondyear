# Matches b stars with GALEX photons

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky
import matplotlib.gridspec as gridspec

##############################################################
##############################################################
# TO DO
# 2/10
# - Stack brightest stars, look for dx dy offset
##############################################################
##############################################################

##############################################################
# Load B stars and make NUV cut for 10-100 cts/s (for 20-15 mag)
##############################################################
bstar = fits.open('../bstar.fits')[1].data
bcut1 = np.where((bstar.nuv_cps > 10.) & (bstar.nuv_cps < 100.))
bstar = bstar[bcut1]
bgal = SkyCoord(bstar.ra*u.degree, bstar.dec*u.degree, frame='icrs').galactic

##############################################################
# Load CSV files to match up times (cut to one file for now), fix gl values
##############################################################
path = "../corrcsv/"
lines = np.loadtxt("../filelists/csv_galcorr.txt", dtype = 'string')

field = fits.open('../corrcsv/00005_0001galcorr.csv.fits')[1].data

#for skyfield in range(len(lines)):
#	field = fits.open(path+lines[skyfield])

##############################################################
# Cut into one second slices & limit bstar data range
##############################################################
names,times,nphotons,stargl,stargb,phogl,phogb = [],[],[],[],[],[],[]
timestep = np.arange(np.round(field.time[-1])+1)
tstep = 210
maxtime = 0

#fig = plt.figure(figsize=(8,8))

#grid = gridspec.GridSpec(5,10,wspace=0.0,hspace=0.0)
gridcount = 0
expsize = 5.
arcminlim = 1.
totdx,totdy = [],[]

# Limit to the first three timesteps to make it run faster
while tstep < timestep[-1]:
	print tstep
	tstep += expsize

	fieldlim = np.where((field.time >= tstep) & (field.time < tstep+expsize))
	fgl = field.gl[fieldlim]
	fgb = field.gb[fieldlim]

	# Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
	for i in range(len(fgl)):
		if min(fgl) < 10. and fgl[i] > 310.:
			fgl[i] = fgl[i] - 360.
	
	# Convert to galactic coordinates and limit the bstar range 
	fgal = SkyCoord(fgl.tolist()*u.degree, fgb.tolist()*u.degree, frame='galactic')
	bcut2 = np.where((bgal.l.degree > np.min(fgl)) & (bgal.l.degree < np.max(fgl)) & (bgal.b.degree > -10.) & (bgal.b.degree < 10.))
	bstar2 = bstar[bcut2]
	bgal2 = SkyCoord(bstar2.ra*u.degree, bstar2.dec*u.degree, frame='icrs').galactic

	##############################################################
	# Find photons within 2' of each star 
	##############################################################	
	starind, photind, starsep, stardist = search_around_sky(bgal2,fgal,arcminlim*u.arcmin)
	bstar2 = bstar2[np.unique(starind)]
	starname = bstar2.galex_id
	bgal2 = SkyCoord(bstar2.ra*u.degree, bstar2.dec*u.degree, frame='icrs').galactic
	bgl = bgal2.l.degree
	bgb = bgal2.b.degree

	fgl = fgal.l.degree
	fgb = fgal.b.degree
	
	# Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
	for i in range(len(fgl)):
		if min(fgl) < 10. and fgl[i] > 310.:
			fgl[i] = fgl[i] - 360.
	
	# Make a list for each star/photon list pair...probably more efficient way to do this
	starnumber = []
	for i in range(len(starind)-1):
		if starind[i] != starind[i+1]:
			starnumber.append(i) 		# Do this so we can match indices with photind

	print starnumber

	pllist = []
	pblist = []
	pllist.append(fgl[photind[0:starnumber[0]+1]])
	pblist.append(fgb[photind[0:starnumber[0]+1]])
	for i in range(len(starnumber)-1):
		pllist.append(fgl[photind[starnumber[i]+1:starnumber[i+1]+1]])
		pblist.append(fgb[photind[starnumber[i]+1:starnumber[i+1]+1]])
	pllist.append(fgl[photind[starnumber[-1]+1:-1]])
	pblist.append(fgb[photind[starnumber[-1]+1:-1]])

	maxnuv = np.argsort(bstar2.nuv_cps[-10:])[::-1] 
	print maxnuv

	##############################################################
	# Get offset for each point 
	##############################################################
	# First translate each star coords to (0,0)
	dx,dy,sumdx,sumdy = [],[],[],[]
	for i in maxnuv:
		dx.append(bgl[i] - pllist[i])
		dy.append(bgb[i] - pblist[i])

	for i in range(len(dx)):
		sumdx = np.concatenate((sumdx,dx[i].tolist()))
		sumdy = np.concatenate((sumdy,dy[i].tolist()))
	
	totdx.append(np.mean(sumdx))
	totdy.append(np.mean(sumdy))

	##############################################################
	# Plot histograms!
	##############################################################

	fig = plt.figure(figsize=(8,8))
	plt.hist(sumdx,bins=50)
	plt.savefig('tstep-'+str(tstep)+'_dx_mean-'+str(np.mean(sumdx))[:8]+'.png')
	plt.close()
	
	fig = plt.figure(figsize=(8,8))
	plt.hist(sumdy,bins=50)
	plt.savefig('tstep-'+str(tstep)+'_dy_mean-'+str(np.mean(sumdy))[:8]+'.png')
	plt.close()

	'''
	plotgrid = 0
	if plotgrid == 1:
		for i in range(len(maxnuv)):
			ax = plt.Subplot(fig,grid[gridcount])
			ax.scatter(bgl[maxnuv[i]],bgb[maxnuv[i]],s=100,facecolor='none',edgecolor='blue',	linewidth='2')
			ax.scatter(pllist[maxnuv[i]],pblist[maxnuv[i]],marker='.',s=1)
			ax.set_xticks([])
			ax.set_yticks([])
			axislim = 0.1 
			ax.text(bgl[maxnuv[i]]-axislim+0.01,bgb[maxnuv[i]]-axislim,str(starname[maxnuv[i]]))
			ax.set_xlim(bgl[maxnuv[i]]-axislim,bgl[maxnuv[i]]+axislim)
			ax.set_ylim(bgb[maxnuv[i]]-axislim,bgb[maxnuv[i]]+axislim)
			fig.add_subplot(ax)
			gridcount+=1
	'''

	otherstuff = 0
	if otherstuff == 1:
		'''
		bskyfix = 0
		if bskyfix == 1: # Do the same gl fix as for the fgl values
			for i in range(len(bgl)):
				if min(bgl) < 10. and bgl[i] > 310.:
					bgl[i] = bgl[i] - 360.
		
		# Check to see that you're not totally screwing up by plotting
		plot = 0
		if plot == 1:
			plt.scatter(bgl,bgb,c='green',s=20,facecolor='none')
			plt.scatter(fgl[photind],fgb[photind],c='orange',s=10,edgecolor='none')
			#plt.scatter(fgl,fgb,c='red',s=3,edgecolor='none')
			plt.show()
			break
	
		# Add photon times for each index
		photon = []
		for i in starind:
			photon.append(sum(starind == i))
	
		#Add values to arrays so you can then make a big table
		cutnames = np.unique(starname)
		bgl = np.unique(bgl)
		bgb = np.unique(bgb)
		fieldgl = fgl[photind].tolist()
		fieldgb = fgb[photind].tolist()
	
		names = np.concatenate((names,starname),axis=0)
		stargl = np.concatenate((stargl,bgl),axis=0)
		stargb = np.concatenate((stargb,bgb),axis=0)
		phogl.append(fieldgl)
		phogb.append(fieldgb)
		times = np.concatenate((times,np.ones(len(names))*tstep),axis=0)
		nphotons = np.concatenate((nphotons,photon),axis=0)
		'''


print totdx
print ' '
print totdy

fig = plt.figure(figsize=(8,8))
plt.scatter(np.arange(0,timestep[-1],5),totdx)
plt.savefig('dxvstime.png')


fig = plt.figure(figsize=(8,8))
plt.scatter(np.arange(0,timestep[-1],5),totdy)
plt.savefig('dyvstime.png')



a1 = np.array(np.arange(0,timestep[-1],5))
a2 = np.array(totdx)
a3 = np.array(totdy)
cols = []
col1 = fits.Column(name='time', format='E', array=a1)
col2 = fits.Column(name='dx', format='E', array=a2)
col3 = fits.Column(name='dy', format='E', array=a3)
cols.append(col1)
cols.append(col2)
cols.append(col3)
new_cols = fits.ColDefs(cols)
hdu = fits.BinTableHDU.from_columns(new_cols)
hdu.writeto('staroffset.fits')


'''
a = []
b = []
for i in range(len(dx)):
	for j in range(len(dx[i])):
		if (np.sqrt(dx[i][j]**2 + dy[i][j]**2) < max(np.abs(dx[i]))/2.):
			a.append(j)
	b.append(a)
	a = []


for i in range(len(dx)):
	plt.scatter(dx[i],dy[i],s=1,alpha=0.5)

#plt.title('Stacked photons 0-50s')


#plt.xlabel('10 Brightest stars, left to right')
#plt.ylabel('time, top to bottom')

#plt.subplots_adjust(hspace=None,wspace=None)



plt.show()
#plt.savefig('bstar_lim_'+str(arcminlim)+'_exp_'+str(expsize)+'.png')
'''
##############################################################
# Now find average photon displacement from each star
##############################################################


##############################################################
# Make the table
##############################################################
make_table = 0
if make_table == 1:
	a1 = np.array(names)
	a2 = np.array(stargl)
	a3 = np.array(stargb)
	a4 = np.array(phogl)
	a5 = np.array(phogb)
	a6 = np.array(times)
	a7 = np.array(nphotons)
	cols = []
	col1 = fits.Column(name='name', format='E', array=a1)
	col2 = fits.Column(name='star_gl', format='E', array=a2)
	col3 = fits.Column(name='star_gb', format='E', array=a3)
	col4 = fits.Column(name='photon_gl', format='20E', array=a4)
	col5 = fits.Column(name='photon_gb', format='E', array=a5)
	col6 = fits.Column(name='times', format='E', array=a6)
	col7 = fits.Column(name='nphotons', format='E', array=a7)
	cols.append(col1)
	cols.append(col2)
	cols.append(col3)
	cols.append(col4)
	cols.append(col5)
	cols.append(col6)
	cols.append(col7)
	new_cols = fits.ColDefs(cols)
	hdu = fits.BinTableHDU.from_columns(new_cols)
	hdu.writeto('00005_0001galcorr_bstarphotons.fits')


