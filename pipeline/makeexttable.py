from astropy.table import Table
from astropy.io import ascii
tablelim = len(t2)
x = np.ones(tablelim)*10.0
obj = np.arange(tablelim)
for i in range(0,len(t2)+1,20000):
    t = Table([obj[i:i+20000-1],t2.ra[i:i+20000-1],t2.dec[i:i+20000-1],x[i:i+20000-1]],names=('obj','ra','dec','size'))
    ascii.write(t,'t2bstars_lowmag'+str(i/20000)+'.txt')


ii = 0
ext = np.loadtxt('extinction_'+str(ii)+'_lowmag.txt',skiprows=13,comments='|',unpack=False)

ra = []
dec = []
cutout_size = []
E_B_V_SandF = []
mean_E_B_V_SandF = []
stdev_E_B_V_SandF = []
max_E_B_V_SandF = []
min_E_B_V_SandF = []
AV_SandF = []
E_B_V_SFD = []
mean_E_B_V_SFD = []
stdev_E_B_V_SFD = []
max_E_B_V_SFD = []
min_E_B_V_SFD = []
AV_SFD = []

for i in range(len(ext)):
    ra.append(ext[i][0])
    dec.append(ext[i][1])
    cutout_size.append(ext[i][2])
    E_B_V_SandF.append(ext[i][3])
    mean_E_B_V_SandF.append(ext[i][4])
    stdev_E_B_V_SandF.append(ext[i][5])
    max_E_B_V_SandF.append(ext[i][6])
    min_E_B_V_SandF.append(ext[i][7])
    AV_SandF.append(ext[i][8])
    E_B_V_SFD.append(ext[i][9])
    mean_E_B_V_SFD.append(ext[i][10])
    stdev_E_B_V_SFD.append(ext[i][11])
    max_E_B_V_SFD.append(ext[i][12])
    min_E_B_V_SFD.append(ext[i][13])
    AV_SFD.append(ext[i][14])

extgal = SkyCoord(ra*u.degree,dec*u.degree, frame='icrs').galactic
gl = extgal.l.degree
gb = extgal.b.degree

a1 = np.array(ra)
a2 = np.array(dec)
a3 = np.array(cutout_size)
a4 = np.array(E_B_V_SandF)
a5 = np.array(mean_E_B_V_SandF)
a6 = np.array(stdev_E_B_V_SandF)
a7 = np.array(max_E_B_V_SandF)
a8 = np.array(min_E_B_V_SandF)
a9 = np.array(AV_SandF)
a10 =np.array(E_B_V_SFD)
a11 =np.array(mean_E_B_V_SFD)
a12 =np.array(stdev_E_B_V_SFD)
a13 =np.array(max_E_B_V_SFD)
a14 =np.array(min_E_B_V_SFD)
a15 =np.array(AV_SFD)
a16 = gl
a17 = gb

col1 = fits.Column(name='ra', format='E', array=a1)
col2 = fits.Column(name='dec', format='E', array=a2)
col3 = fits.Column(name='cutout_size', format='E', array=a3)
col4 = fits.Column(name='E_B_V_SandF', format='E', array=a4)
col5 = fits.Column(name='mean_E_B_V_SandF', format='E', array=a5)
col6 = fits.Column(name='stdev_E_B_V_SandF', format='E', array=a6)
col7 = fits.Column(name='max_E_B_V_SandF', format='E', array=a7)
col8 = fits.Column(name='min_E_B_V_SandF', format='E', array=a8)
col9 = fits.Column(name='AV_SandF', format='E', array=a9)
col10 = fits.Column(name='E_B_V_SFD', format='E', array=a10)
col11 = fits.Column(name='mean_E_B_V_SFD', format='E', array=a11)
col12 = fits.Column(name='stdev_E_B_V_SFD', format='E', array=a12)
col14 = fits.Column(name='max_E_B_V_SFD', format='E', array=a14)
col13 = fits.Column(name='min_E_B_V_SFD', format='E', array=a13)
col15 = fits.Column(name='AV_SFD', format='E', array=a15)
col16 = fits.Column(name='gl', format='E', array=a16)
col17 = fits.Column(name='gb', format='E', array=a17)

new_cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17])
hdu = fits.BinTableHDU.from_columns(new_cols)
hdu.writeto('extinction'+str(ii)+'_lowmag.fits')



