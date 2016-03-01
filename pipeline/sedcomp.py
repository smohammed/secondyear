import numpy as np
import pysynphot
from astropy.io import fits, ascii
from astropy.table import Table, hstack, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import wcs
import os
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

pickles = Table.read('picklemags_laphare_final.txt', format='ascii')
pa = fits.open('pickles_paper.fits')[1].data

nuvmag, Bmag, Vmag, Umag, Rmag, Imag, Jmag, Hmag, Kmag, umag, gmag, rmag, imag = [], [], [], [], [], [], [], [], [], [], [], [], []
for ind in range(1, 132):
    star_file = os.path.join(os.environ['PYSYN_CDBS'], 'grid/pickles/dat_uvi/', 'pickles_'+str(ind)+'.fits')
    star = pysynphot.FileSpectrum(star_file)
    #obs_wfc3 = pysynphot.Observation(star, pysynphot.ObsBandpass('wfc3,uvis1,f814w'), force='extrap')
    obs_star_nuv = pysynphot.Observation(star, pysynphot.ObsBandpass('galex,nuv'))

    obs_star_U = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,u'))
    obs_star_B = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,b'))
    obs_star_V = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,v'))
    obs_star_R = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,r'))
    obs_star_I = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,i'), force='extrap')

    #obs_star_J = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,j'))
    #obs_star_H = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,h'))
    #obs_star_K = pysynphot.Observation(star, pysynphot.ObsBandpass('johnson,k'))

    obs_star_u = pysynphot.Observation(star, pysynphot.ObsBandpass('sdss,u'))
    obs_star_g = pysynphot.Observation(star, pysynphot.ObsBandpass('sdss,g'))
    obs_star_r = pysynphot.Observation(star, pysynphot.ObsBandpass('sdss,r'))
    obs_star_i = pysynphot.Observation(star, pysynphot.ObsBandpass('sdss,i'))

    nuvmag.append(obs_star_nuv.effstim('ABmag'))

    Umag.append(obs_star_U.effstim('ABmag'))
    Bmag.append(obs_star_B.effstim('ABmag'))
    Vmag.append(obs_star_V.effstim('ABmag'))
    Rmag.append(obs_star_R.effstim('ABmag'))

    #Jmag.append(obs_star_J.effstim('ABmag'))
    #Hmag.append(obs_star_H.effstim('ABmag'))
    #Kmag.append(obs_star_K.effstim('ABmag'))

    umag.append(obs_star_u.effstim('ABmag'))
    gmag.append(obs_star_g.effstim('ABmag'))
    rmag.append(obs_star_r.effstim('ABmag'))
    imag.append(obs_star_i.effstim('ABmag'))

#binset=obs_wfc3.binwave)

nuvmag = np.array(nuvmag)

Umag = np.array(Umag)
Bmag = np.array(Bmag)
Vmag = np.array(Vmag)
Rmag = np.array(Rmag)

Jmag = np.array(Jmag)
Hmag = np.array(Hmag)
Kmag = np.array(Kmag)

umag = np.array(umag)
gmag = np.array(gmag)
rmag = np.array(rmag)
imag = np.array(imag)


obs_star_U.effstim('ABmag') - obs_star_B.effstim('ABmag')

