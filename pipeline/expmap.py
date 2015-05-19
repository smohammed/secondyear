from astropy.io import fits
import numpy as np

#################################################
# Open files to mask
#################################################
files1 = np.loadtxt('../filelists/aspcorr_new_list.txt', dtype='string')
#files1 = np.loadtxt('../filelists/truncaspcorr_new_list.txt',dtype='string') # For 0-20deg
path = "../scst/"
photonlist = np.loadtxt('../filelists/photontimes.txt', dtype='string')

def expmap(begin,end):
    #################################################
    # Make array exposure map
    #################################################
    # These can be arbitrary
    nx=81.*10
    ny=81.*10

    # add a buffer around edge just for good measure
    #expim=np.zeros((360.0*60+nx,20.0*60+ny))
    expim=np.zeros((200.0*60+nx,200.0*60+ny))

    #################################################
    # Creates a circular mask array of 1s bordered by 0s
    #################################################
    mask = np.zeros((nx,ny))
    radius = (1.24/2)*60.*10

    for i in np.arange(nx):
        for j in np.arange(ny):
            a = i-400.
            b = j-400.
            if np.sqrt(a**2+b**2) < radius:
                mask[i][j] = 1.

    #################################################
    # Files too big to make one map, make mini maps instead. First select files
    #################################################
    startgl = 0
    endgl = 0
    while (float(files1[startgl].split('_')[3])/10.) < begin:
        startgl += 1
    while (float(files1[endgl].split('_')[3])/10.) < end:
        endgl += 1

    for line in range(startgl, endgl, 1):                         # Iterate through scst files
        print line
        print float(files1[line].split('_')[3])/10.
        d = fits.open(path+files1[line])[1].data                    # Open exp file

        for pholine in range(len(photonlist)):                      # Find correct photon times file
            if files1[line].split('_')[3] == photonlist[pholine][:5]:
                pfile = pholine

        photimes = np.loadtxt('../times/'+photonlist[pfile],comments='#',delimiter=',',unpack=False,dtype=  str)[:-1]                           # Load correct photon times file

        delpho = []
        for i in range(len(photimes)-1):
           delpho.append(np.where((d['T']-d['T'][0])==float(photimes[i])))

        q = []
        for i in range(len(delpho)):
            if len(delpho[i][0]) > 0:
                q.append(delpho[i][0][0])

        #################################################
        # Convert RA, DEC to gl, gb then to pixels
        #################################################
        ##### CONVERT RA, DEC to GL, GB, and then to PIXELS, assuming that each pixel is 1 arcmin (thus     the x60 below)  
        ##### THESE ARE ARRAYS of SCAN position at a given time (in 1 second intervals)plt.
        #gl = SkyCoord(d.ra_acs*u.degree, d.dec_acs*u.degree, frame='icrs').galactic.l.degree
        #gb = SkyCoord(d.ra_acs*u.degree, d.dec_acs*u.degree, frame='icrs').galactic.b.degree
    
        gl = d.SCAN_GL[q]
        gb = d.SCAN_GB[q]
    
        # Convert to arcmin then to pixels, 12000/(20*60) pixels/arcmin
        gx = (gl * 60.) * 12000/(20*60)
        gy = ((gb + 10.) * 60.) * 12000/(20*60)
        ''' 
        for i in range(len(gx)):
            if gx[i] < 0:
                gx[i] = 0
            elif gx[i] > 360*60.:
                gx[i] = 360*60.

        for i in range(len(gy)):
            if gy[i] < 0:
                gy[i] = 0
            elif gy[i] > 20*60.:
                gy[i] = 20*60.
        '''
        for i in range(len(gx)):
            if gx[i] < 0:
                gx[i] = 0
            elif gx[i] > 200*60.:
                gx[i] = 200*60.

        for i in range(len(gy)):
            if gy[i] < 0:
                gy[i] = 0
            elif gy[i] > 200*60.:
                gy[i] = 200*60.

        ##### ADD BORDER TO PIXEL VALUE
        # use 40 offset to keep from trailing off to account for errors at the border of the map
        gx = gx + 400.
        gy = gy + 400.

        #################################################
        # Apply mask and only take data with > 15000 photon counts
        #################################################

        ##### ONLY USE SCAN POSITIONS FROM TIMES WHERE WE RECEIVED MORE THAN 15000 photons per second.      You can do something similar by only taking scan positions that have
        ##### photons that you used to build the images.

        for j in range(len(d[q])-1):
            if d.NDCTEC[q][j] > 15000.:
                if np.shape(expim[gx[j]-405.:gx[j]+405., gy[j]-405.:gy[j]+405.]) == np.shape(mask):
                    expim[gx[j]-405.:gx[j]+405., gy[j]-405.:gy[j]+405.] = expim[gx[j]-405.:gx[j]+405., gy[j]-405.:gy[j]+405.] + mask
                else:
                    print 'shape of expim != shape of mask'

    hdu = fits.PrimaryHDU(expim)
    print 'finished'
    return hdu.writeto('../expmap_12000_'+str(float(files1[startgl].split('_')[3])/10.)+'to'+str(   float(files1[endgl-1].split('_')[3])/10.)+'.fits')

for glstep in range(0, 100, 20):
    expmap(glstep,glstep+20)
