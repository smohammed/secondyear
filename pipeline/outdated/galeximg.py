from astropy.io import fits
import numpy as np
from numpy import loadtxt

# For Dun
# It's been a while since I used this (probably 2015) so a lot of this may be confusing
# For the inputs, I just gather all the files within that range (for a 20x20 deg block)
# Including the up and down runs for each scan and near the end (line 195) I use
# np.histogram2d with all gl, gb photon counts and bin it using linspace for each side
# Then I transpose it so the image comes out in the correct orientation and then finally
# I save it as a fits file. The binnum was 12000 and this was too large for my computer.
# You can experiment with the bin number but I feel like you'll have a better idea than me.

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
#a = 10             # Start gl
#b = 20             # End gl
#c = -10            # Start gb
#d = 10             # End gb


def intensitymap(a, b, c, d, pointing, correction):
    ############################################
    # Set the file pathways
    ############################################

    if pointing == 1:
        path = "../intmapfiles/"
        lines = loadtxt("../filelists/pointingtimes.txt", comments="#", delimiter=",", unpack=False,   dtype=str)

    if pointing == 2:
        path = "../corrcsv/"
        #lines=loadtxt("../filelists/AIS_GAL_SCAN_csv.txt", comments="#", delimiter=",", unpack=False, dtype=str)
        lines = loadtxt("../filelists/csv_galcorr.txt", dtype='string')

    if correction == 1:
        offpath = "../dxdyoffsets/"
        offlines = loadtxt("../filelists/dxdycorrectionfiles.txt", dtype='string')

    ############################################
    # Grab appropriate files using gl range
    ############################################
    startgl = 0
    endgl = 0

    if pointing == 1:
        while (float(lines[startgl].split('_')[0])/10.) < a:
            startgl += 1
        while (float(lines[endgl].split('_')[0])/10.) < b:
            endgl += 1
        endgl = endgl - 1

        photons = fits.open(path+lines[startgl])[1].data
        photons2 = fits.open(path+lines[startgl+1])[1].data

        gl = np.concatenate((photons.gl_cor, photons2.gl_cor))
        gb = np.concatenate((photons.gb_cor, photons2.gb_cor))

    if pointing == 2:
        while (float(lines[startgl].split('_')[0])/10.) < a:
            startgl += 1
        while (float(lines[endgl].split('_')[0])/10.) < b:
            endgl += 1
        endgl = endgl - 1

        photons = fits.open(path+lines[startgl])[1].data
        photlim1 = np.where((photons.gb >= c) & (photons.gb < d))
        gl = photons.gl[photlim1]
        gb = photons.gb[photlim1]
        time = photons.time[photlim1]

        gl[(min(gl) < 10.) & (gl > 350.)] = gl[(min(gl) < 10.) & (gl > 350.)] - 360.

        photlim2 = np.where((gl >= a) & (gl < b))

        gl = gl[photlim2]
        gb = gb[photlim2]
        time = time[photlim2]

        if correction == 1:
            offset = fits.open(offpath+offlines[startgl])[1].data

            for i in range(len(offset)-1):
                offsetind = np.where((time >= offset.time[i]) & (time < offset.time[i+1]))
                gl[offsetind] = gl[offsetind] - offset.dx[i]/60.
                gb[offsetind] = gb[offsetind] - offset.dy[i]/60.

    print 'startgl = ', startgl
    print 'endgl  = ', endgl

    ############################################
    # Combine photons from up and down runs, stack data horizontally
    ############################################
    if pointing == 1:
        for i in range(startgl+2, endgl, 2):
            print i
            try:
                f1 = fits.open(path+lines[i])[1].data
                f1gl_cor = f1.gl_cor
                f1gb_cor = f1.gb_cor

            except IOError:                             # if the file is corrupted, skip it
                print 'Index '+str(i)+' is corrupted'

            try:
                f2 = fits.open(path+lines[i+1])[1].data
                f2gl_cor = f2.gl_cor
                f2gb_cor = f2.gb_cor

            except IOError:                             # if the file is corrupted, skip it
                print 'Index '+str(i+1)+' is corrupted'

            glinit = np.concatenate((f1gl_cor, f2gl_cor))
            gbinit = np.concatenate((f1gb_cor, f2gb_cor))
            gl = np.hstack([gl, glinit])
            gb = np.hstack([gb, gbinit])

    if pointing == 2:
        for i in range(startgl+1, endgl+1):
            print str(i) + ' for file ' + lines[i]
            try:
                f1 = fits.open(path+lines[i])[1].data
                f1lim1 = np.where((f1.gb >= c) & (f1.gb < d))
                f1gl = f1.gl[f1lim1]
                f1gb = f1.gb[f1lim1]
                time = f1.time[f1lim1]
                f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] = f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] - 360.

                f1lim2 = np.where((f1gl >= a) & (f1gl < b))
                f1gl = f1gl[f1lim2]
                f1gb = f1gb[f1lim2]
                time = time[f1lim2]

                if correction == 1:
                    offset = fits.open(offpath+offlines[i])[1].data
                    for i in range(len(offset)-1):
                        offsetind = np.where((time >= offset.time[i]) & (time < offset.time[i+1]))
                        f1gl[offsetind] = f1gl[offsetind] - offset.dx[i]/60.
                        f1gb[offsetind] = f1gb[offsetind] - offset.dy[i]/60.

            except IOError:                             # if the file is corrupted, skip it
                print 'Index '+str(i)+' is corrupted'
            except ValueError:
                print 'No data available for ' + lines[i]
            gl = np.hstack([gl, f1gl])
            gb = np.hstack([gb, f1gb])

    # Be sure to write the last scan!
    if pointing == 2:
        if float(lines[endgl+1].split('_')[0])/10. == 359.6:
            print 'Writing last slice 359.6 deg'
            b = b + 1

            f1 = fits.open(path+lines[endgl + 1])[1].data
            photlim1 = np.where((photons.gb >= c) & (photons.gb < d))
            f1gl = f1.gl[f1lim1]
            f1gb = f1.gb[f1lim1]
            time = f1.time[f1lim1]

            f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] = f1gl[(min(f1gl) < 10.) & (f1gl > 350.)] - 360.
            f1lim2 = np.where((f1gl >= a) & (f1gl < b))

            f1gl = f1gl[f1lim2]
            f1gb = f1gb[f1lim2]
            time = time[f1lim2]

            gl = np.hstack([gl, f1gl])
            gb = np.hstack([gb, f1gb])

    ############################################
    # Bin data and plot stuff
    ############################################
    if pointing == 1:
        binnum = 1200

        H, xbins, ybins = np.histogram2d(gl, gb, bins=(np.linspace(a, b, binnum), np.linspace(c, d,   binnum)))
        #fig=plt.figure(figsize=(10,10))
        #ax=plt.axes()

        fits.PrimaryHDU(H.T).writeto('../intensitymap'+str(binnum)+'_'+str(a)+'_'+str(b)+'.fits')

    if pointing == 2:
        binnum = 12000
        print binnum
        H, xbins, ybins = np.histogram2d(gl, gb, bins=(np.linspace(a, b, binnum), np.linspace(c, d, binnum)))

        H = H.T
        '''
        # characterize your data in terms of a linear translation from XY pixels to gl, gb

        # lambda function given min, max, n_pixels, return spacing, middle value.
        linwcs=lambda x, y, n: ((x-y)/n, (x+y)/2)

        cdeltaX, crvalX=linwcs(np.amin(gl), np.amax(gl), len(gl))
        cdeltaY, crvalY=linwcs(np.amin(gb), np.amax(gb), len(gb))

        # wcs code ripped from
        # http://docs.astropy.org/en/latest/wcs/index.html

        w=wcs.WCS(naxis=2)

        # what is the center pixel of the XY grid.
        w.wcs.crpix=[len(gl)/2, len(gb)/2]

        # what is the galactic coordinate of that pixel.
        w.wcs.crval=[crvalX, crvalY]

        # what is the pixel scale in lon, lat.
        w.wcs.cdelt=np.array([cdeltaX, cdeltaY])

        # you would have to determine if this is in fact a tangential projection.
        w.wcs.ctype=["GLON", "GLAT"]

        # write the HDU object WITH THE HEADER
        header=w.to_header()
        '''

        if correction == 0:
            #return fits.PrimaryHDU(H,header=header).writeto('../intmapcsvcorr'+str(binnum)+'_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'.fits')
            return H, gl, gb, xbins, bins
        elif correction == 1:
            return fits.PrimaryHDU(H, header=header).writeto('../intmapcsvcorr'+str(binnum)+'_gl_'+str(a)+'to'+str(b)+'_gb_'+str(c)+'to'+str(d)+'_corr.fits')

