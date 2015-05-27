from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

def plotmap(startgl, endgl, startgb, endgb, bins):
    bins = bins
    a = startgl     # Start gl
    b = endgl       # End gl

    c = startgb
    d = endgb

    # Load files
    if bins == 1200:
        intlist = np.loadtxt('../filelists/intmapcsvcorr1200.txt', dtype='string')
        explist = np.loadtxt('../filelists/expmaps1200.txt', dtype='string')

        for i in range(len(intlist)):
            if intlist[i] == 'intmapcsvcorr1200_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits':
                fileindex = i

        ###################################################################
        # Set a scale so we can cut the correct portion of the sky from the big expmap 
        ###################################################################
        pixscalex = (b - 20.)/20.
        pixscaley = (d)
        offset = 0 #1 # Found for bins=1200 via testing

        ###################################################################
        # Load the 20x20 photon and exposure maps
        ###################################################################
        phomap = np.sqrt(fits.open('../photonmaps1200/intmapcsvcorr1200_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits')[0].data)
        expmap = fits.open('../expmaps1200/'+explist[fileindex])[0].data.T

        ###################################################################
        # Set pixel limits on the image to remove the edges of the expmap, then apply to the map
        ###################################################################
        ymin = 39  # is gb = -10
        ymax = 1238  # is gb = 10
        xmin = 39  # +1199 * pixscalex + offset
        xmax = 1238  # +1199 * pixscalex + offset

    if bins == 12000:
        intlist = np.loadtxt('../filelists/intmapcsvcorr12000.txt', dtype='string')
        explist = np.loadtxt('../filelists/expmaps12000.txt', dtype='string')

        for i in range(len(intlist)):
            if intlist[i] == 'intmapcsvcorr12000_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits':
                fileindex = i

        ###################################################################
        # Load the 20x20 deg photon and exposure maps
        ###################################################################
        phomap = np.sqrt(fits.open('../photonmaps12000/intmapcsvcorr12000_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits')[0].data)
        expmap = fits.open('../expmaps12000/'+explist[fileindex])[0].data.T

        ###################################################################
        # Set pixel limits on the image to remove the edges of the expmap, then apply to the map
        ###################################################################
        ymin = 405  # is gb = -10
        ymax = 12810-406  # is gb = 10
        xmin = 405  # +1199 * pixscalex + offset
        xmax = 12810-406  # +1199 * pixscalex + offset

    #expmap1 = expmap[39:1238,39+1199*pixscalex:1238+1199*pixscalex]
    #expmap1 = expmap[39:1238,39+offset+1199*pixscalex:1238+offset+1199*pixscalex]
    expmap1 = expmap[ymin:ymax, xmin:xmax]

    ###################################################################
    # Combine the exposure map with the intensity map
    ###################################################################
    intmap = phomap/expmap1
    intmap[np.isnan(intmap)] = 0
    intmap[intmap == np.inf] = 0
    intmap[intmap == -1*np.inf] = 0

    # Also plot bstars?
    bstarplot = 0
    if bstarplot == 1:
        ###################################################################
        # Find gl,gb of bstars in that region, restrict by flux
        ###################################################################
        bstars = fits.open('../bstar.fits')[1].data
        gl = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.l.degree
        gb = SkyCoord(bstars.ra*u.degree, bstars.dec*u.degree, frame='icrs').galactic.b.degree
        blim = np.where((bstars.nuv_cps > cutfluxlow) & (bstars.nuv_cps < cutfluxhigh) & (gl > a) & (gl < b) & (gb > -10) & (gb < 10))
        gl = gl[blim]
        gb = gb[blim]

        ###################################################################
        # Now plot!
        ###################################################################
        #plt.scatter(gl,gb,facecolor='none',edgecolor='red')
        #plt.xlabel('gl')
        #plt.ylabel('gb')

    if bins == 1200:
        #return plt.imshow(intmap,vmin=0,vmax=0.7,origin='lower',extent=[a,b,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray), plt.show()
        hdu = fits.PrimaryHDU(intmap)
        return hdu.writeto('../comb1200_gl'+str(a)+'to'+str(b)+'.fits')

    if bins == 12000:
        #return plt.imshow(intmap,vmin=0,vmax=0.05,origin='lower',extent=[a,b,-10,10],interpolation='nearest',aspect='auto',cmap=cm.gray), plt.show()
        hdu = fits.PrimaryHDU(intmap)
        return hdu.writeto('../comb12000_gl'+str(a)+'to'+str(b)+'.fits')

'''
for i in range(0, 140, 20):
    plotmap(i, i+20, -10, 10, 12000)
    print i

for i in range(160, 240, 20):
    plotmap(i, i+20, -10, 10, 12000)
    print i
'''
for i in range(280, 360, 20):
    print 'gl = '+ str(i)
    plotmap(i, i+20, -10, 10, 12000)


#plotmap(340, 360, -10, 10, 12000)