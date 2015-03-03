# Matches b stars with GALEX photons

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky
import matplotlib.gridspec as gridspec

def plotdata(x,y):
	from matplotlib.ticker import NullFormatter

	# the random data
	nullfmt   = NullFormatter()         # no labels
	
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02
	
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	
	# start with a rectangular Figure
	plt.figure(1, figsize=(8,8))
	
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)
	
	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	axHisty.yaxis.set_major_formatter(nullfmt)
	
	# the scatter plot:
	axScatter.scatter(x, y)
	
	# now determine nice limits by hand:
	binwidth = 0.25
	xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
	lim = ( int(xymax/binwidth) + 1) * binwidth
	
	axScatter.set_xlim( (-lim, lim) )
	axScatter.set_ylim( (-lim, lim) )
	
	bins = np.arange(-lim, lim + binwidth, binwidth)
	axHistx.hist(x, bins=bins)
	axHisty.hist(y, bins=bins, orientation='horizontal')
	
	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	
	return plt.show()

##############################################################
##############################################################
# TO DO
# 2/19
# - Reset dx,dy on star peak, take some radius from that, compute mean, std
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

skyfield = 150 #what skyfield?

field = fits.open('../corrcsv/'+lines[skyfield])[1].data
print lines[skyfield]

#for skyfield in range(len(lines)):
#	field = fits.open(path+lines[skyfield])

##############################################################
# Cut into one second slices & limit bstar data range
##############################################################
names,times,nphotons,stargl,stargb,phogl,phogb = [],[],[],[],[],[],[]
timestep = np.arange(np.round(field.time[-1])+1)
tstep = 0
expsize = 5.
arcminlim = 2.
totdx,totdy,totphot2, totphot15,stdsumdx,stdsumdy = [],[],[],[],[],[]
timesteplimit = -1 	#Index of the timestep limits. -1 is the entire set
tempdgl, tempdgb = 0.,0.
fix = 0.

print 'timesteps = ', timestep[timesteplimit]

plotgrid = 0
if plotgrid == 1:
	fig = plt.figure(figsize=(8,8))
	grid = gridspec.GridSpec(5,10,wspace=0.0,hspace=0.0)
	gridcount = 0


while (tstep < timestep[timesteplimit-5]):# and (tstep != timestep[timesteplimit-5]):
	print 'tstep = ', tstep
	tstep += expsize

	fieldlim = np.where((field.time >= tstep) & (field.time < tstep+expsize))
	fgl = field.gl[fieldlim]
	fgb = field.gb[fieldlim]

	# Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
	if len(fgl) > 0.:
		fgl[(min(fgl) < 10.) & (fgl > 350.)] = fgl[(min(fgl) < 10.) & (fgl > 350.)] - 360.

		# Convert to galactic coordinates and limit the bstar range 
		fgal = SkyCoord(fgl.tolist()*u.degree, fgb.tolist()*u.degree, frame='galactic')
		print 'len(bgal)= ', len(bgal)
		print 'len(fgal)= ', len(fgl)
	
		bcut2 = np.where((bgal.l.degree > np.min(fgl)) & (bgal.l.degree < np.max(fgl)) & (bgal.b.degree > -10.) & (bgal.b.degree < 10.))
		bstar2 = bstar[bcut2]
		bgal2 = SkyCoord(bstar2.ra*u.degree, bstar2.dec*u.degree, frame='icrs').galactic
	
		##############################################################
		# Find photons within 2' of each star 
		##############################################################	
		starind, photind, starsep, stardist = search_around_sky(bgal2,fgal,arcminlim*u.arcmin)
		bstar2am = bstar2[np.unique(starind)]
		starname = bstar2am.galex_id
		bgal2am = SkyCoord(bstar2am.ra*u.degree, bstar2am.dec*u.degree, frame='icrs').galactic
		bgl = bgal2am.l.degree
		bgb = bgal2am.b.degree
		fgl = fgal.l.degree
		fgb = fgal.b.degree
		
		# Fix neg values of fgl again so that are ~360 deg to neg numbers to keep circular shape
		fgl[(min(fgl) < 10.) & (fgl > 350.)] = fgl[(min(fgl) < 10.) & (fgl > 350.)] - 360.
	
		# Make a list for each star/photon list pair...probably more efficient way to do this
		starnumber = []
		for i in range(len(starind)-1):
			if starind[i] != starind[i+1]:
				starnumber.append(i) 		# Do this so we can match indices with photind	
	
		if len(starnumber) > 2.:
			pllist = []
			pblist = []
			pllist.append(fgl[photind[0:starnumber[0]+1]]-tempdgl)
			pblist.append(fgb[photind[0:starnumber[0]+1]]-tempdgb)
			for i in range(len(starnumber)-1):
				pllist.append(fgl[photind[starnumber[i]+1:starnumber[i+1]+1]])
				pblist.append(fgb[photind[starnumber[i]+1:starnumber[i+1]+1]])
			pllist.append(fgl[photind[starnumber[-1]+1:-1]])
			pblist.append(fgb[photind[starnumber[-1]+1:-1]])
		
			maxnuv = np.argsort(bstar2am.nuv_cps[-10:])[::-1]
			#print 'maxnuv =', maxnuv
		
			##############################################################
			# Get offset for each point, recenter, calc mean, stdev
			##############################################################
			# First translate each star coords to (0,0)
			dx,dy,sumdx,sumdy = [],[],[],[]
			for i in maxnuv:
				dx.append((bgl[i] - pllist[i])*60.) # dx and dy in arcminutes
				dy.append((bgb[i] - pblist[i])*60.)
		
			# Sum photons from all stars 
			for i in range(len(dx)):
				sumdx = np.concatenate((sumdx,dx[i].tolist()))
				sumdy = np.concatenate((sumdy,dy[i].tolist()))
			
			bns = 50 
	
			tx = np.histogram(sumdx,bins=bns)[1][np.argmax(np.histogram(sumdx,bins=bns)[0])]
			ty = np.histogram(sumdy,bins=bns)[1][np.argmax(np.histogram(sumdy,bins=bns)[0])]
	
			sumdx = sumdx - tx
			sumdy = sumdy - ty
	
			#plotdata(sumdx,sumdy)
			print tx
			print ty
	
			radlim = np.where(np.sqrt(sumdx**2 + sumdy**2) < 1.)
			
			sumdx1am = sumdx[radlim]
			sumdy1am = sumdy[radlim]
			
			totdx.append(np.mean(sumdx1am) + tempdgl)
			totdy.append(np.mean(sumdy1am) + tempdgb)
			
			if fix == 1.:
				tempdgl = np.mean(sumdx1am)
				tempdgb = np.mean(sumdy1am)
			stdsumdx.append(np.std(sumdx1am)/np.sqrt(len(sumdx1am)))
			stdsumdy.append(np.std(sumdy1am)/np.sqrt(len(sumdy1am)))
	
			# Add up all photons within 30 arcsec
			radlim30as = np.where(np.sqrt(sumdx**2 + sumdy**2) < 1/2.)
			totphot2.append(len(sumdx[radlim30as]))
	
			print tempdgl
			print tempdgb
			##############################################################
			# Make plots
			##############################################################
			plot_hist = 0.
			if plot_hist == 1:
				fig = plt.figure(figsize=(8,8))
				plt.hist(sumdx,bins=50)
				plt.savefig('tstep-'+str(tstep)+'_dx_mean-'+str(np.mean(sumdx))[:8]+'.png')
				plt.close()
		
				fig = plt.figure(figsize=(8,8))
				plt.hist(sumdy,bins=50)
				plt.savefig('tstep-'+str(tstep)+'_dy_mean-'+str(np.mean(sumdy))[:8]+'.png')
				plt.close()
	
			if plotgrid == 1:
				for i in range(len(maxnuv)):
					ax = plt.Subplot(fig,grid[gridcount])
					ax.scatter(bgl[maxnuv[i]],bgb[maxnuv[i]],s=100,facecolor='none',edgecolor='blue',		linewidth='2')
					ax.scatter(pllist[maxnuv[i]],pblist[maxnuv[i]],marker='.',s=1)
					ax.set_xticks([])
					ax.set_yticks([])
					axislim = 0.1 
					ax.text(bgl[maxnuv[i]]-axislim+0.01,bgb[maxnuv[i]]-axislim,str(starname[maxnuv[i]])	)
					ax.set_xlim(bgl[maxnuv[i]]-axislim,bgl[maxnuv[i]]+axislim)
					ax.set_ylim(bgb[maxnuv[i]]-axislim,bgb[maxnuv[i]]+axislim)
					fig.add_subplot(ax)
					gridcount+=1
		
		else:
			print 'no stars in slice tstep =', tstep
			totdx.append(0.)
			totdy.append(0.)
			stdsumdx.append(0.)
			stdsumdy.append(0.)
			totphot2.append(0.)
	
	
		find_Nphot15 = 0
		if find_Nphot15 == 1:
			##############################################################
			# Find photons in a 15' radius now
			##############################################################
			starind15, photind15, starsep15, stardist15 = search_around_sky(bgal2,fgal,15.*u.arcmin)
			bstar215 = bstar2[np.unique(starind15)]
			starname15 = bstar215.galex_id
			bgal215 = SkyCoord(bstar215.ra*u.degree, bstar215.dec*u.degree, frame='icrs').galactic
			bgl15 = bgal215.l.degree
			bgb15 = bgal215.b.degree
		
			fgl15 = fgal.l.degree
			fgb15 = fgal.b.degree
			
			# Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
			fgl15[(min(fgl15) < 10.) & (fgl15 > 350.)] = fgl15[(min(fgl15) < 10.) & (fgl15 > 350.)] - 	360.
		
			# Make a list for each star/photon list pair...probably more efficient way to do this
			starnumber15 = []
			for i in range(len(starind15)-1):
				if starind15[i] != starind15[i+1]:
					starnumber15.append(i) 		# Do this so we can match indices with photind
		
			#print 'starnumber =', starnumber
		
			if len(starnumber15) > 0.:
				pllist15 = []
				pblist15 = []
				pllist15.append(fgl15[photind15[0:starnumber15[0]+1]])
				pblist15.append(fgb15[photind15[0:starnumber15[0]+1]])
				for i in range(len(starnumber15)-1):
					pllist15.append(fgl15[photind15[starnumber15[i]+1:starnumber15[i+1]+1]])
					pblist15.append(fgb15[photind15[starnumber15[i]+1:starnumber15[i+1]+1]])
				pllist15.append(fgl15[photind15[starnumber15[-1]+1:-1]])
				pblist15.append(fgb15[photind15[starnumber15[-1]+1:-1]])
			
				maxnuv15 = np.argsort(bstar215.nuv_cps[-10:])[::-1] 
				print 'maxnuv15 =', maxnuv15
			
				##############################################################
				# Get offset for each point, recenter, calc mean, stdev
				##############################################################
				# First translate each star coords to (0,0)
				dx15,sumdx15 = [],[]
				for i in maxnuv:
					dx15.append((bgl15[i] - pllist15[i])*60.)
			
				for i in range(len(dx15)):
					sumdx15 = np.concatenate((sumdx15,dx15[i].tolist()))
				
				totphot15.append(len(sumdx15))

	else:
		print 'len(fgl) = 0'
		totdx.append(0.)
		totdy.append(0.)
		stdsumdx.append(0.)
		stdsumdy.append(0.)
		totphot2.append(0.)

totdx.append(0.)
totdy.append(0.)
stdsumdx.append(0.)
stdsumdy.append(0.)
totphot2.append(0.)

print 'totdx =', len(totdx)
print 'totdy =', len(totdy)
print 'totphot2 =', len(totphot2)

f, axarr = plt.subplots(3, sharex=True,figsize=(16,10))
axarr[0].errorbar(np.arange(0,timestep[timesteplimit],5)-1,totdx,yerr=stdsumdx,fmt='o')
axarr[0].set_ylabel('dx [arcmin]')
axarr[0].set_xlim(0,timestep[timesteplimit])
axarr[1].errorbar(np.arange(0,timestep[timesteplimit],5)-1,totdy,yerr=stdsumdy,fmt='o')
axarr[1].set_ylabel('dy [arcmin]')
axarr[2].scatter(np.arange(0,timestep[timesteplimit],5)-1, totphot2)
axarr[2].set_xlabel('Time [s]')
axarr[2].set_ylabel('NPhot 2arcmin')
axarr[2].axhline(y=100)
f.subplots_adjust(hspace=0)
if fix == 1.:
	f.savefig('../dxdyplots/dxdyvst'+str(lines[skyfield])[:5]+'_fix.png')
else:
	f.savefig('../dxdyplots/dxdyvst'+str(lines[skyfield])[:5]+'_nofix.png')
#plt.show()

##############################################################
# Make a table of dx, dy, Nphotons in 2' and 15' 
##############################################################
make_star_table = 0
if make_star_table == 1.:
	a1 = np.array(np.arange(0,timestep[timesteplimit],5))
	a2 = np.array(totdx)
	a3 = np.array(totdy)
	a4 = np.array(totphot2)
	a5 = np.array(totphot15)
	cols = []
	col1 = fits.Column(name='time', format='E', array=a1)
	col2 = fits.Column(name='dx', format='E', array=a2)
	col3 = fits.Column(name='dy', format='E', array=a3)
	col4 = fits.Column(name='Nphot2', format='E', array=a4)
	col5 = fits.Column(name='Nphot15', format='E', array=a5)	
	cols.append(col1)
	cols.append(col2)
	cols.append(col3)
	cols.append(col4)
	cols.append(col5)
	new_cols = fits.ColDefs(cols)
	hdu = fits.BinTableHDU.from_columns(new_cols)
	hdu.writeto('staroffset.fits')

##############################################################
# Now find average photon displacement from each star
##############################################################


