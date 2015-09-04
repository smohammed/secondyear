import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['figure.figsize'] = 16, 8
matplotlib.rcParams['font.size'] = 20

nuvjvsjk = 1
nuvbvsbv = 0
jhvshk = 0
grvsnuvg = 0


#star = Table.read('newfield_gal_2mass_t2_jlim_13.5_3arcsec.txt', format='ascii')

star = Table.read('newfield_gal_ipac_2mass_matches_3arcsec_j_lt13.5_tests.txt',format='ascii')
newt = Table.read('galex0data_2mass_t2.txt', format='ascii')
pickles = Table.read('picklemags_laphare_test.txt', format='ascii')

# Set range
scut = np.where((star['gb_sex'] > -10) & (star['gb_sex'] < -5))
scut2 = np.where((star['gb_sex'] > -5) & (star['gb_sex'] < 0))
scut3 = np.where((star['gb_sex'] > 0) & (star['gb_sex'] < 5))
scut4 = np.where((star['gb_sex'] > 5) & (star['gb_sex'] < 10))
ncut = np.where((newt['gb_galex'] > -10) & (newt['gb_galex'] < -5))
ncut2 = np.where((newt['gb_galex'] > -5) & (newt['gb_galex'] < 0))
ncut3 = np.where((newt['gb_galex'] > 0) & (newt['gb_galex'] < 5))
ncut4 = np.where((newt['gb_galex'] > 5) & (newt['gb_galex'] < 10))


# GALEX, Pickles, SExtractor order
nuv = ['nuv_mag', 'nuv', 'nuv', 2.9720]
b = ['BJmag', 'b', 'BJmag', 1.3429]
v = ['VJmag', 'v', 'VJmag', 1.0]
j = ['j', 'j', 'j', 0.2876]
h = ['h', 'h', 'h', 0.1783]
k = ['k', 'k', 'k', 0.1170]
u = ['u','u_AB','u_AB',1.5916]
g = ['g','g_AB','g_AB',1.1838]
r = ['r','r_AB','r_AB',0.8664]
ib  = ['i','i_AB','i_AB',0.6418]

if nuvjvsjk == 1:
    x1 = j
    x2 = k
    y1 = nuv
    y2 = j
    extx = -0.1
    exty = 8
    #extx = 0.08
    #exty = 3.5
    extdx = x1[3] - x2[3]
    extdy = y1[3] - y2[3]
if nuvbvsbv == 1:
    x1 = b
    x2 = v
    y1 = nuv
    y2 = b
    extx = 1.4
    exty = 3
    #extx = 0.2
    #exty = 3.4
    extdx = x1[3] - x2[3]
    extdy = y1[3] - y2[3]
if jhvshk == 1:
    x1 = h
    x2 = k
    y1 = j
    y2 = h
    extx = 0.25
    exty = -0.1
    extdx = x1[3] - x2[3]
    extdy = y1[3] - y2[3]

if grvsnuvg == 1:
    x1 = g
    x2 = r
    y1 = nuv
    y2 = g
    extx = 4
    exty = -0.175
    extdx = x1[3] - x2[3]
    extdy = y1[3] - y2[3]


# Open Figure
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
plt.rc('legend', **{'fontsize': 15})

# SExtractor
a1 = ax1.scatter(star[x1[2]][scut]-star[x2[2]][scut], star[y1[2]][scut]-star[y2[2]][scut], edgecolor='none', alpha=0.3)
b1 = ax2.scatter(star[x1[2]][scut2]-star[x2[2]][scut2], star[y1[2]][scut2]-star[y2[2]][scut2], edgecolor='none', alpha=0.3)
c1 = ax3.scatter(star[x1[2]][scut3]-star[x2[2]][scut3], star[y1[2]][scut3]-star[y2[2]][scut3], edgecolor='none', alpha=0.3)
d1 = ax4.scatter(star[x1[2]][scut4]-star[x2[2]][scut4], star[y1[2]][scut4]-star[y2[2]][scut4], edgecolor='none', alpha=0.3)


# GALEX data
a2 = ax1.scatter(newt[x1[0]][ncut]-newt[x2[0]][ncut], newt[y1[0]][ncut]-newt[y2[0]][ncut], facecolor='red', edgecolor='none', s=30, alpha=0.3)
b2 = ax2.scatter(newt[x1[0]][ncut2]-newt[x2[0]][ncut2], newt[y1[0]][ncut2]-newt[y2[0]][ncut2], facecolor='red', edgecolor='none', s=30, alpha=0.3)
c2 = ax3.scatter(newt[x1[0]][ncut3]-newt[x2[0]][ncut3], newt[y1[0]][ncut3]-newt[y2[0]][ncut3], facecolor='red', edgecolor='none', s=30, alpha=0.3)
d2 = ax4.scatter(newt[x1[0]][ncut4]-newt[x2[0]][ncut4], newt[y1[0]][ncut4]-newt[y2[0]][ncut4], facecolor='red', edgecolor='none', s=30, alpha=0.3)

#Pickles
a3 = ax1.scatter(pickles[x1[1]]-pickles[x2[1]], pickles[y1[1]]-pickles[y2[1]], c='black', s=10)
b3 = ax2.scatter(pickles[x1[1]]-pickles[x2[1]], pickles[y1[1]]-pickles[y2[1]], c='black', s=10)
c3 = ax3.scatter(pickles[x1[1]]-pickles[x2[1]], pickles[y1[1]]-pickles[y2[1]], c='black', s=10)
d3 = ax4.scatter(pickles[x1[1]]-pickles[x2[1]], pickles[y1[1]]-pickles[y2[1]], c='black', s=10)


for j in range(0, len(pickles), 10):
    ax1.annotate(pickles['name'][j], xy=((pickles[x1[1]][j]-pickles[x2[1]][j])+0.01, (pickles[y1[1]][j]-pickles[y2[1]][j])+0.07), size=15)
    ax2.annotate(pickles['name'][j], xy=((pickles[x1[1]][j]-pickles[x2[1]][j])+0.01, (pickles[y1[1]][j]-pickles[y2[1]][j])+0.07), size=15)
    ax3.annotate(pickles['name'][j], xy=((pickles[x1[1]][j]-pickles[x2[1]][j])+0.01, (pickles[y1[1]][j]-pickles[y2[1]][j])+0.07), size=15)
    ax4.annotate(pickles['name'][j], xy=((pickles[x1[1]][j]-pickles[x2[1]][j])+0.01, (pickles[y1[1]][j]-pickles[y2[1]][j])+0.07), size=15)




# Add Extinction vector
ax1.arrow(extx, exty, extdx, extdy, head_length=0.1,head_width=0.07,color='black')
ax2.arrow(extx, exty, extdx, extdy, head_length=0.1,head_width=0.07,color='black')
ax3.arrow(extx, exty, extdx, extdy, head_length=0.1,head_width=0.07,color='black')
ax4.arrow(extx, exty, extdx, extdy, head_length=0.1,head_width=0.07,color='black')

# Set all labels
ax1.set_title('3", J < 13.5, -10 < gb < -5')
#ax1.set_title('6", -10 < gb < -5')
ax2.set_title('-5 < gb < 0')
ax3.set_title('0 < gb < 5')
ax4.set_title('5 < gb < 10')

if nuvjvsjk == 1:
    ax1.set_xlim((-0.25, 2.5))
    ax1.set_ylim((-2, 16))
    ax2.set_xlim((-0.25, 2.5))
    ax2.set_ylim((-2, 16))
    ax3.set_xlim((-0.25, 2.5))
    ax3.set_ylim((-2, 16))
    ax4.set_xlim((-0.25, 2.5))
    ax4.set_ylim((-2, 16))
    ax1.set_ylabel('NUV - J')
    ax3.set_xlabel('J - K')
    ax3.set_ylabel('NUV - J')
    ax4.set_xlabel('J - K')

if nuvbvsbv == 1:
    ax1.set_xlim((-0.5, 2))
    ax1.set_ylim((2, 10))
    ax2.set_xlim((-0.5, 2))
    ax2.set_ylim((2, 10))
    ax3.set_xlim((-0.5, 2))
    ax3.set_ylim((2, 10))
    ax4.set_xlim((-0.5, 2))
    ax4.set_ylim((2, 10))
    ax1.set_ylabel('NUV - B')
    ax3.set_xlabel('B - V')
    ax3.set_ylabel('NUV - B')
    ax4.set_xlabel('B - V')

if jhvshk == 1:
    ax1.set_xlim((-0.1,0.4))
    ax1.set_ylim((-0.2,1.0))
    ax2.set_xlim((-0.1,0.4))
    ax2.set_ylim((-0.2,1.0))
    ax3.set_xlim((-0.1,0.4))
    ax3.set_ylim((-0.2,1.0))
    ax4.set_xlim((-0.1,0.4))
    ax4.set_ylim((-0.2,1.0))
    ax1.set_ylabel('J - H')
    ax3.set_xlabel('H - K')
    ax3.set_ylabel('J - H')
    ax4.set_xlabel('H - K')

if grvsnuvg == 1:
    ax1.set_xlim((-3,8))
    ax1.set_ylim((-1.,2.0))
    ax2.set_xlim((-3,8))
    ax2.set_ylim((-1.,2.0))
    ax3.set_xlim((-3,8))
    ax3.set_ylim((-1.,2.0))
    ax4.set_xlim((-3,8))
    ax4.set_ylim((-1.,2.0))
    ax1.set_ylabel('NUV - g (AB mag)')
    ax3.set_xlabel('g - r (AB mag)')
    ax3.set_ylabel('NUV - g (AB mag)')
    ax4.set_xlabel('g - r (AB mag)')

if nuvjvsjk == 1:
    ax1.legend([a1, a2, a3], ['SExtractor', 'GALEX', 'Pickles'], scatterpoints=1, loc=4)
    #ax1.legend([a1, a2], ['SExtractor','GALEX'], scatterpoints=1, loc=1)

if nuvbvsbv == 1:
    ax1.legend([a1, a2, a3], ['SExtractor', 'GALEX', 'Pickles'], scatterpoints=1, loc=2)  
    #ax1.legend([a1, a2], ['SExtractor','GALEX'], scatterpoints=1, loc=2)

if jhvshk == 1:
    ax1.legend([a1, a2, a3], ['SExtractor', 'GALEX', 'Pickles'], scatterpoints=1, loc=2)

if grvsnuvg == 1:
    ax1.legend([a1, a2, a3], ['SExtractor', 'GALEX', 'Pickles'], scatterpoints=1, loc=4)


# Reference line
if nuvjvsjk == 1:
    ax1.axhline(y=4,color='black')
    ax2.axhline(y=4,color='black')
    ax3.axhline(y=4,color='black')
    ax4.axhline(y=4,color='black')

f.subplots_adjust(wspace=0)
plt.show()



