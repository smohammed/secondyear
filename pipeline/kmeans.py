import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from sklearn.cluster import KMeans
from sklearn import preprocessing
matplotlib.rcParams['figure.figsize'] = 10,10


def kmeans(file,t0,t1,kpts):
    field = fits.open('../corrcsv/00005_0001galcorr.csv.fits')[1].data
    fieldlim = np.where((field.time >=t0) & (field.time < t1))
    gl = field.gl[fieldlim]
    gb = field.gb[fieldlim]
    
    for i in range(len(gl)):
        if min(gl) < 10. and gl[i] > 310.:
            gl[i] = gl[i] - 360.
    
    binnum = 1200
    grid = np.vstack([gl,gb]).T
    H, glbins, gbbins = np.histogram2d(gl,gb,bins=binnum)
    
    scaler = preprocessing.StandardScaler()
    n_clusters = kpts
    clf = KMeans(n_clusters)
    clf.fit(scaler.fit_transform(grid))
    cluster_centers = scaler.inverse_transform(clf.cluster_centers_)
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot()
    
    # plot density
    ax = plt.axes()
    ax.imshow(H.T, origin='lower',vmin=0.0,vmax=6., interpolation='nearest', aspect='auto',extent=[glbins[0], glbins[-1],gbbins[0], gbbins[-1]])
    
    # plot cluster centers
    cluster_centers = scaler.inverse_transform(clf.cluster_centers_)
    ax.scatter(cluster_centers[:, 0], cluster_centers[:, 1],s=40, c='w', edgecolors='k')
    
    # plot cluster boundaries
    gl_centers = 0.5 * (glbins[1:] + glbins[:-1])
    gb_centers = 0.5 * (gbbins[1:] + gbbins[:-1])
    
    Xgrid = np.meshgrid(gl_centers, gb_centers)
    Xgrid = np.array(Xgrid).reshape((2, binnum * binnum)).T
    
    H = clf.predict(scaler.transform(Xgrid)).reshape((binnum, binnum))
    
    for i in range(n_clusters):
        Hcp = H.copy()
        flag = (Hcp == i)
        Hcp[flag] = 1
        Hcp[~flag] = 0
    
        ax.contour(gl_centers, gb_centers, Hcp, [-0.5, 0.5],linewidths=1, colors='k')
    
    #ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
    #ax.set_xlim(glbins[0], glbins[-1])
    #ax.set_ylim(gbbins[0], gbbins[-1])
    ax.set_xlim(-0.5,1.5)
    ax.set_ylim(8.5,11.)
    
    ax.set_xlabel('gl')
    ax.set_ylabel('gb')
    
    return plt.show() 

t0 = 0.0
k = 500

kmeans(0,t0,t0+1,k)
#kmeans(0,t0+1,t0+2,k)
#kmeans(0,t0+2,t0+3,k)