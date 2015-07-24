from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

# Compare sextractor number counts and catalog match counts with various plots

f1counts = Table.read('../newfield.sex_total_05-68.txt',format='ascii')
f1matches = Table.read('../starcatalog_05-68_2mass_t2.txt',format='ascii')

print 'f1counts = ', len(f1counts)
print 'f1matches = ', len(f1matches)



h1 = plt.hist(f1counts['nuv'],bins=100,color='red',label='counts')
h2 = plt.hist(f1matches['nuv'],bins=100,color='blue',label='matches')
plt.title('New field Sextractor counts vs catalog matches')
plt.legend(loc=2)

