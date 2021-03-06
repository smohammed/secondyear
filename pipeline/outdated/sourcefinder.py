# Matches b stars with GALEX photons

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matching import search_around_sky
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from plotdata import plotdata
import os


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

##############################################################
##############################################################
# TO DO
# stuff
##############################################################
##############################################################

##############################################################
# Load B stars and make NUV cut for 10-100 cts/s (for 20-15 mag)
##############################################################
#bstar = fits.open('../bstar.fits')[1].data
bstar = fits.open('../starcat.fits')[1].data
bcut1 = np.where((bstar.nuv_cps > 5.) & (bstar.nuv_cps < 100.))
bstar = bstar[bcut1]
bgal = SkyCoord(bstar.gl*u.degree, bstar.gb*u.degree, frame='galactic')

##############################################################
# Load CSV files to match up times (cut to one file for now), fix gl values
##############################################################
path = "../corrcsv/"
lines = np.loadtxt("../filelists/csv_galcorr.txt", dtype='string')

fix = 0.
#skyfield = 11  # what skyfield?

# 0,20,60,78,100,120,130,150,180,190,210,220,240,280,300,330,360,390,420,440

#field = fits.open(path+lines[skyfield])[1].data
#print lines[skyfield]

for skyfield in range(353, len(lines), 1):
    print lines[skyfield]
    field = fits.open(path+lines[skyfield])[1].data

    ##############################################################
    # Cut into one second slices & limit bstar data range
    ##############################################################
    timestep = np.arange(np.round(field.time[-1])+1)
    timesteplimit = -1  # Index of the timestep limits. -1 is the entire set

    plotgrid = 0
    if plotgrid == 1:
        fig = plt.figure(figsize=(8, 8))
        grid = gridspec.GridSpec(5, 10, wspace=0.0, hspace=0.0)
        gridcount = 0

    # Apply offset from last iteration?
    off = 0.
    for offsetrun in range(1):
        tstep = 0
        expsize = 5
        arcminlim = 2.
        totdx, totdy, totphot2,  totphot15, stdsumdx, stdsumdy = [], [], [], [], [], []
        timesteplimit = -1  # Index of the timestep limits. -1 is the entire set
        tempdgl, tempdgb = 0., 0.
        timelist = np.arange(0, timestep[timesteplimit-5], 5)

        # Find regions where the tx, ty offset moves the data far from 0,0
        print 'timesteps = ', timestep[timesteplimit]
        dxdyerror = []

        if off == 1.:
            offsetpath = '../dxdyoffsets/'
            if os.path.exists('../dxdyoffsets/'+str(lines[skyfield])[:10]+'offset.fits'):
                if offsetrun == 0:
                    os.remove('../dxdyoffsets/'+str(lines[skyfield])[:10]+'offset.fits')
                offsetfile = fits.open(offsetpath+lines[skyfield][:10]+'offset.fits')[1].data
            else:
                pass
        print 'offsetrun =', offsetrun
        
        while (tstep < timestep[timesteplimit-5]):
            print ' '
            print 'tstep = ', tstep
    
            fieldlim = np.where((field.time >= tstep) & (field.time < tstep+expsize))
            fgl = field.gl[fieldlim]
            fgb = field.gb[fieldlim]
    
            if len(fgl) > 0.:  # Removes areas with no photon data
                # Fix neg values of fgl that are ~360 deg to neg numbers to keep circular shape
                fgl[(min(fgl) < 10.) & (fgl > 350.)] = fgl[(min(fgl) < 10.) & (fgl > 350.)] - 360.
    
                # Convert to galactic coordinates and limit the bstar range
                fgal = SkyCoord(fgl.tolist()*u.degree, fgb.tolist()*u.degree, frame='galactic')
    
                bcut2 = np.where((bgal.l.degree > np.min(fgl)) & (bgal.l.degree < np.max(fgl)) & (bgal.b.degree > -10.) & (bgal.b.degree < 10.))
                bstar2 = bstar[bcut2]
                bgal2 = SkyCoord(bstar2.gl*u.degree, bstar2.gb*u.degree, frame='galactic')
                print len(bgal2)
                ##############################################################
                # Find photons within 2' of each star
                ##############################################################
                starind, photind, starsep, stardist = search_around_sky(bgal2, fgal, arcminlim*u.arcmin)
                starnumber = []
                if len(starind) > 0.:
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
                    for i in range(len(starind)-1):
                        if starind[i] != starind[i+1]:
                            starnumber.append(i)        # Do this so we can match indices with photind
    
                if len(starnumber) > 2.:
                    pllist = []
                    pblist = []
                    pllist.append(fgl[photind[0:starnumber[0]+1]])
                    pblist.append(fgb[photind[0:starnumber[0]+1]])
                    for i in range(len(starnumber)-1):
                        pllist.append(fgl[photind[starnumber[i]+1:starnumber[i+1]+1]])
                        pblist.append(fgb[photind[starnumber[i]+1:starnumber[i+1]+1]])
                    pllist.append(fgl[photind[starnumber[-1]+1:-1]])
                    pblist.append(fgb[photind[starnumber[-1]+1:-1]])
    
                    maxnuv = np.argsort(bstar2am.nuv_cps[-10:])[::-1]
    
                    ##############################################################
                    # Get offset for each point, recenter, calc mean, stdev
                    ##############################################################
                    # First translate each star coords to (0,0)
                    dx, dy, sumdx, sumdy = [], [], [], []
                    for i in maxnuv:
                        dx.append((bgl[i] - pllist[i])*60.)  # dx and dy in arcminutes
                        dy.append((bgb[i] - pblist[i])*60.)
    
                    # Sum photons from all stars
                    for i in range(len(dx)):
                        sumdx = np.concatenate((sumdx, dx[i].tolist()))
                        sumdy = np.concatenate((sumdy, dy[i].tolist()))
    
                    # Original star positions
                    if ((tstep % 50 == 0.) & (offsetrun != 0.)):
                        plotdata(sumdx, sumdy)
    
                    '''
                    # Apply offset from last iteration, then add it to the photons
                    if ((off == 1.) & (offsetrun != 0.)):
                        offind = np.where(offsetfile.time == (tstep)) # Remember to change to 2.5
                        print 'offind =', offind[0][0]
                        print 'dxoffset = ', offsetfile.dx[offind][0]
                        print 'dyoffset = ', offsetfile.dy[offind][0]
    
                        sumdx = sumdx - offsetfile.dx[offind][0]
                        sumdy = sumdy - offsetfile.dy[offind][0]
                    '''
    
                    # Find peak of a histogram for dx,dy, remove peaks that are too offset from 0,0
                    bns = 25
    
                    # Define some test data which is close to Gaussian
                    histx, bin_edgesx = np.histogram(sumdx, bins=bns, density=True)
                    bin_centersx = (bin_edgesx[:-1] + bin_edgesx[1:])/2
    
                    histy, bin_edgesy = np.histogram(sumdy, bins=bns, density=True)
                    bin_centersy = (bin_edgesy[:-1] + bin_edgesy[1:])/2
    
                    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
                    px0 = [1., 0., 1.]
                    py0 = [1., 0., 1.]
    
                    try:
                        coeffx, var_matrixx = curve_fit(gauss, bin_centersx, histx, p0=px0)
                        coeffy, var_matrixy = curve_fit(gauss, bin_centersy, histy, p0=py0)
                        tx = coeffx[1]
                        ty = coeffy[1]
    
                        '''
                        print 'Fitted x mean = ', coeffx[1]
                        print 'Fitted x standard deviation = ', coeffx[2]
                        print 'Fitted y mean = ', coeffy[1]
                        print 'Fitted y standard deviation = ', coeffy[2]
                        '''
    
                    except RuntimeError:
                        tx = 0.
                        ty = 0.
                        print 'tx = ', tx
                        print 'ty = ', ty
    
                    #tx = np.histogram(sumdx,bins=bns)[1][np.argmax(np.histogram(sumdx,bins=bns)[0])]
                    #ty = np.histogram(sumdy,bins=bns)[1][np.argmax(np.histogram(sumdy,bins=bns)[0])]
    
                    # Positions with correction applied
                    if ((tstep % 50 == 0.) & (offsetrun != 0.)):
                        plotdata(sumdx, sumdy)
    
                    # Make sure tx, ty offset isn't too large. If it is, make a note
                    #if ((tx < 2.) & (ty < 2.)):
                    if (np.sqrt(tx**2 + ty**2) < 2.):
                        print 'tx =', tx
                        print 'ty =', ty
                        pass
                    else:
                        dxdyerror.append(tstep)
                        print 'txty error at', tstep
                        tx,  ty = 0., 0.
    
                    # Recenter dx, dy using histogram
                    if fix == 1:
                        sumdx = sumdx - tx
                        sumdy = sumdy - ty
    
                    # Position after histogram recentering
                    if ((tstep % 50 == 0.) & (offsetrun != 0.)):
                        plotdata(sumdx, sumdy)
    
                    # Remove photons that are outside 1"
                    #plotdata(sumdx,sumdy)
                    radlim = np.where(np.sqrt(sumdx**2 + sumdy**2) < 1.)
    
                    sumdx1am = sumdx[radlim]
                    sumdy1am = sumdy[radlim]
    
                    #print 'tempdgl =',tempdgl
                    #print 'mean sumdx1am =', np.mean(sumdx1am)
    
                    # temp coords apply a fix from the last position to the next one
                    totdx.append(np.mean(sumdx1am))  # + tempdgl)
                    totdy.append(np.mean(sumdy1am))  # + tempdgb)
    
                    tempdgl = np.mean(sumdx1am)
                    tempdgb = np.mean(sumdy1am)
    
                    stdsumdx.append(np.std(sumdx1am)/np.sqrt(len(sumdx1am)))
                    stdsumdy.append(np.std(sumdy1am)/np.sqrt(len(sumdy1am)))
    
                    # Add up all photons within 30 arcsec
                    radlim30as = np.where(np.sqrt(sumdx**2 + sumdy**2) < 1/2.)
                    radlim2am = np.where(np.sqrt(sumdx**2 + sumdy**2) < 2.)
                    totphot2.append(len(sumdx[radlim2am]))
    
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
                           ax.scatter(bgl[maxnuv[i]],bgb[maxnuv[i]],s=100,facecolor='none',edgecolor='blue',       linewidth='2')
                           ax.scatter(pllist[maxnuv[i]],pblist[maxnuv[i]],marker='.',s=1)
                           ax.set_xticks([])
                           ax.set_yticks([])
                           axislim = 0.1 
                           ax.text(bgl[maxnuv[i]]-axislim+0.01,bgb[maxnuv[i]]-axislim,str(starname[maxnuv[i]]) )
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
                    fgl15[(min(fgl15) < 10.) & (fgl15 > 350.)] = fgl15[(min(fgl15) < 10.) & (fgl15 > 350.)] -   360.
                
                    # Make a list for each star/photon list pair...probably more efficient way to do this
                    starnumber15 = []
                    for i in range(len(starind15)-1):
                        if starind15[i] != starind15[i+1]:
                            starnumber15.append(i)      # Do this so we can match indices with photind
                
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
    
            else:  # If there are no photons in that region
                print 'len(fgl) =', len(fgl)
                totdx.append(0.)
                totdy.append(0.)
                stdsumdx.append(0.)
                stdsumdy.append(0.)
                totphot2.append(0.)
            tstep += expsize
    
        totdx.append(0.)
        totdy.append(0.)
        stdsumdx.append(0.)
        stdsumdy.append(0.)
        totphot2.append(0.)
    
        print 'totdx = ', len(totdx)
        print 'totdy = ', len(totdy)
    
        print 'totphot2 =', len(totphot2)
        print 'dxdyerror =', dxdyerror
        print len(timelist)
        timelist = np.append(timelist, timestep[timesteplimit])
    

        f, axarr = plt.subplots(3, sharex=True, figsize=(16, 10))
        axarr[0].errorbar(timelist, totdx, yerr=stdsumdx, fmt='o')
        axarr[0].set_ylabel('dx [arcmin]')
        axarr[0].set_xlim(0, timestep[timesteplimit])
        axarr[0].set_title(str(lines[skyfield][:10])+', tstep = 0.')
        axarr[1].errorbar(timelist, totdy, yerr=stdsumdy, fmt='o')
        axarr[1].set_ylabel('dy [arcmin]')
        axarr[2].scatter(timelist,  totphot2)
        axarr[2].set_xlabel('Time [s]')
        axarr[2].set_ylabel('NPhot 2arcmin')
        axarr[2].axhline(y=100)
        f.subplots_adjust(hspace=0)
    
        if off == 0:
            if fix == 1.:
                f.savefig('../dxdyplots/obsnuv/dxdyvst'+str(lines[skyfield])[:10]+'obsnuv_fix.png')
            else:
                f.savefig('../dxdyplots/obsnuv/dxdyvst'+str(lines[skyfield])[:10]+'obsnuv_nofix.png')
            #plt.show()
        elif off == 1:
            if fix == 1.:
                f.savefig('../dxdyplots/dxdyvst'+str(lines[skyfield])[:10]+'_fix_iter10.png')
            else:
                f.savefig('../dxdyplots/dxdyvst'+str(lines[skyfield])[:10]+'_nofix_iter10.png')
            #plt.show()
        
        ##############################################################
        # Make a table of dx, dy, Nphotons in 2' and 15'
        ##############################################################
        make_star_table = 1
        if make_star_table == 1.:
            a1 = np.array(np.arange(0, timestep[timesteplimit], 5))
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
    
            if off == 0.:
                hdu.writeto('../dxdyoffsets/obsnuv/'+str(lines[skyfield])[:10]+'offset_obsnuv.fits')
    
            if (off == 1.):  # ((off == 1.) & (offsetrun == 0.)):
                if os.path.exists('../dxdyoffsets/obsnuv/'+str(lines[skyfield])[:10]+'offset_obsnuv.fits'):
                    os.remove('../dxdyoffsets/obsnuv/'+str(lines[skyfield])[:10]+'offset_obsnuv.fits')
    
                hdu.writeto('../dxdyoffsets/obsnuv/'+str(lines[skyfield])[:10]+'offset_obsnuv.fits')
    