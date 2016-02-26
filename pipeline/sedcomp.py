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



'''
# Filter files
directory = '../PICKLES/filt/'
nuvfilt = Table.read(directory+'galex/NUV.pb', format='ascii')
Bfilt = Table.read(directory+'wfi/B.pb', format='ascii', data_start=1)
Vfilt = Table.read(directory+'wfi/V.pb', format='ascii', data_start=1)
Ufilt = Table.read(directory+'wfi/U.pb', format='ascii', data_start=1)
Rfilt = Table.read(directory+'wfi/R.pb', format='ascii', data_start=1)
Ifilt = Table.read(directory+'wfi/I.pb', format='ascii', data_start=1)

Jfilt = Table.read(directory+'2mass/J.pb', format='ascii', data_start=1)
Hfilt = Table.read(directory+'2mass/H.pb', format='ascii', data_start=1)
Kfilt = Table.read(directory+'2mass/Ks.pb', format='ascii', data_start=1)

ufilt = Table.read(directory+'sdss/up.pb', format='ascii', data_start=1)
gfilt = Table.read(directory+'sdss/gp.pb', format='ascii', data_start=1)
rfilt = Table.read(directory+'sdss/rp.pb', format='ascii', data_start=1)
ifilt = Table.read(directory+'sdss/ip.pb', format='ascii', data_start=1)

# Pysynphot Pickles files
nuvmag, Bmag, Vmag, Umag, Rmag, Imag, Jmag, Hmag, Kmag, umag, gmag, rmag, imag = [], [], [], [], [], [], [], [], [], [], [], [], []

for ind in range(1, 132):
    star_file = os.path.join(os.environ['PYSYN_CDBS'], 'grid/pickles/dat_uvi/', 'pickles_'+str(ind)+'.fits')
    star = pysynphot.FileSpectrum(star_file)

    nuvstarcut = np.where(np.in1d(star.wave, nuvfilt['col1']))
    Bstarcut = np.where(np.in1d(star.wave, Bfilt['col1']))
    Vstarcut = np.where(np.in1d(star.wave, Vfilt['col1']))
    Ustarcut = np.where(np.in1d(star.wave, Ufilt['col1']))
    Rstarcut = np.where(np.in1d(star.wave, Rfilt['col1']))
    Istarcut = np.where(np.in1d(star.wave, Ifilt['col1']))

    Jstarcut = np.where(np.in1d(star.wave, Jfilt['col1']))
    Hstarcut = np.where(np.in1d(star.wave, Hfilt['col1']))
    Kstarcut = np.where(np.in1d(star.wave, Kfilt['col1']))

    ustarcut = np.where(np.in1d(star.wave, ufilt['col1']))
    gstarcut = np.where(np.in1d(star.wave, gfilt['col1']))
    rstarcut = np.where(np.in1d(star.wave, rfilt['col1']))
    istarcut = np.where(np.in1d(star.wave, ifilt['col1']))

    nuvfiltcut = np.where(np.in1d(nuvfilt['col1'], star.wave))
    Bfiltcut = np.where(np.in1d(Bfilt['col1'], star.wave))
    Vfiltcut = np.where(np.in1d(Vfilt['col1'], star.wave))
    Ufiltcut = np.where(np.in1d(Ufilt['col1'], star.wave))
    Rfiltcut = np.where(np.in1d(Rfilt['col1'], star.wave))
    Ifiltcut = np.where(np.in1d(Ifilt['col1'], star.wave))

    Jfiltcut = np.where(np.in1d(Jfilt['col1'], star.wave))
    Hfiltcut = np.where(np.in1d(Hfilt['col1'], star.wave))
    Kfiltcut = np.where(np.in1d(Kfilt['col1'], star.wave))

    ufiltcut = np.where(np.in1d(ufilt['col1'], star.wave))
    gfiltcut = np.where(np.in1d(gfilt['col1'], star.wave))
    rfiltcut = np.where(np.in1d(rfilt['col1'], star.wave))
    ifiltcut = np.where(np.in1d(ifilt['col1'], star.wave))

    nuvmag.append(-2.5*np.log10(np.sum(star.flux[nuvstarcut]*star.wave[nuvstarcut]*nuvfilt['col2'][nuvfiltcut])/(np.sum(3631*3*10**10/star.wave[nuvstarcut]*nuvfilt['col2'][nuvfiltcut]))))
    Bmag.append(-2.5*np.log10(np.sum(star.flux[Bstarcut]*star.wave[Bstarcut]*Bfilt['col2'][Bfiltcut])/(np.sum(3631*3*10**10/star.wave[Bstarcut]*Bfilt['col2'][Bfiltcut]))))
    Vmag.append(-2.5*np.log10(np.sum(star.flux[Vstarcut]*star.wave[Vstarcut]*Vfilt['col2'][Vfiltcut])/(np.sum(3631*3*10**10/star.wave[Vstarcut]*Vfilt['col2'][Vfiltcut]))))
    Umag.append(-2.5*np.log10(np.sum(star.flux[Ustarcut]*star.wave[Ustarcut]*Ufilt['col2'][Ufiltcut])/(np.sum(3631*3*10**10/star.wave[Ustarcut]*Ufilt['col2'][Ufiltcut]))))
    Rmag.append(-2.5*np.log10(np.sum(star.flux[Rstarcut]*star.wave[Rstarcut]*Rfilt['col2'][Rfiltcut])/(np.sum(3631*3*10**10/star.wave[Rstarcut]*Rfilt['col2'][Rfiltcut]))))
    Imag.append(-2.5*np.log10(np.sum(star.flux[Istarcut]*star.wave[Istarcut]*Ifilt['col2'][Ifiltcut])/(np.sum(3631*3*10**10/star.wave[Istarcut]*Ifilt['col2'][Ifiltcut]))))

    Jmag.append(-2.5*np.log10(np.sum(star.flux[Jstarcut]*star.wave[Jstarcut]*Jfilt['col2'][Jfiltcut])/(np.sum(3631*3*10**10/star.wave[Jstarcut]*Jfilt['col2'][Jfiltcut]))))
    Hmag.append(-2.5*np.log10(np.sum(star.flux[Hstarcut]*star.wave[Hstarcut]*Hfilt['col2'][Hfiltcut])/(np.sum(3631*3*10**10/star.wave[Hstarcut]*Hfilt['col2'][Hfiltcut]))))
    Kmag.append(-2.5*np.log10(np.sum(star.flux[Kstarcut]*star.wave[Kstarcut]*Kfilt['col2'][Kfiltcut])/(np.sum(3631*3*10**10/star.wave[Kstarcut]*Kfilt['col2'][Kfiltcut]))))

    umag.append(-2.5*np.log10(np.sum(star.flux[ustarcut]*star.wave[ustarcut]*ufilt['col2'][ufiltcut])/(np.sum(3631*3*10**10/star.wave[ustarcut]*ufilt['col2'][ufiltcut]))))
    gmag.append(-2.5*np.log10(np.sum(star.flux[gstarcut]*star.wave[gstarcut]*gfilt['col2'][gfiltcut])/(np.sum(3631*3*10**10/star.wave[gstarcut]*gfilt['col2'][gfiltcut]))))
    rmag.append(-2.5*np.log10(np.sum(star.flux[rstarcut]*star.wave[rstarcut]*rfilt['col2'][rfiltcut])/(np.sum(3631*3*10**10/star.wave[rstarcut]*rfilt['col2'][rfiltcut]))))
    imag.append(-2.5*np.log10(np.sum(star.flux[istarcut]*star.wave[istarcut]*ifilt['col2'][ifiltcut])/(np.sum(3631*3*10**10/star.wave[istarcut]*ifilt['col2'][ifiltcut]))))

nuvmag = np.array(nuvmag)
Bmag = np.array(Bmag)
Vmag = np.array(Vmag)
Umag = np.array(Umag)
Rmag = np.array(Rmag)
Imag = np.array(Imag)
Jmag = np.array(Jmag)
Hmag = np.array(Hmag)
Kmag = np.array(Kmag)
umag = np.array(umag)
gmag = np.array(gmag)
rmag = np.array(rmag)
imag = np.array(imag)


data = [nuvmag, Umag, Bmag, Vmag, Rmag, Imag, Jmag, Hmag, Kmag, umag, gmag, rmag, imag]
ascii.write(data, '../pysynphot_pickles.txt', names=['nuv', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'u', 'g', 'r', 'i'])
'''