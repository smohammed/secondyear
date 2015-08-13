def plotdata(x, y):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    # the random data
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, edgecolor='None')

    # now determine nice limits by hand:
    binwidth = 0.08
    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    lim = (int(xymax/binwidth) + 1) * binwidth

    #axScatter.set_xlim( (-lim, lim) )
    #axScatter.set_ylim( (-lim, lim) )

    #axScatter.set_xlim((-0.4, 0.6))
    #axScatter.set_ylim((-0.3, 0.9))

    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x, bins=bins)
    axHisty.hist(y, bins=bins, orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axScatter.axhline(y=0., xmin=-20, xmax=20, color='k')
    axScatter.axvline(x=0., ymin=-20, ymax=20, color='k')

    return plt.show()


def ccplot(x1,x2,y1,y2):
    star = fits.open('../starcatalog.fits')[1].data
    pickles = Table.read('../picklemags.txt',format='ascii')
    lims = np.array([0,  20,  40,  60,  80, 100, 120, 160, 180, 200, 220, 280, 300, 320, 340])
    for i in lims:
        scut = np.where((np.abs(star.gb > 5.)) & (star.gl > i) & (star.gl < i+20))
        scut2 = np.where((np.abs(star.gb < 5.)) & (star.gl > i) & (star.gl < i+20))
        f, (ax1, ax2) = plt.subplots(1, 2)
        ax1.scatter(star[x1][scut]-star[x2][scut],star[y1][scut]-star[y2][scut],edgecolor='none',alpha=0.2,c=star.gb[scut],vmin=0,vmax=5)
        ax1.set_title('gb > 5, gl = '+str(i)+' to '+str(i+20))
        #ax1.set_xlim((-0.5,1.5))
        #ax1.set_ylim((3,12.5))
        ax1.set_xlabel(x1+' - '+x2)
        ax1.set_ylabel(y1+' - '+y2)
        ax1.plot([-0.5,1.5],[3,12])
        ax1.scatter(pickles[x1]-pickles[x2],pickles[y1]-pickles[y2],c='black',s=5)
        for i in range(len(pickles)):
            ax1.annotate(pickles['name'][i][:-4],xy=(pickles[x1][i]-pickles[x2][i],pickles[y1][i]-pickles[y2][i]),size=10)

        ax2.scatter(star[x1][scut2]-star[x2][scut2],star[y1][scut2]-star[y2][scut2],edgecolor='none',alpha=0.2,c=star.gb[scut2],vmin=5,vmax=10)
        ax2.set_title('gb < 5, gl = '+str(i)+' to '+str(i+20))
        #ax2.set_xlim((-0.5,1.5))
        #ax2.set_ylim((3,12.5))
        ax2.set_xlabel(x1+' - '+x2)
        ax2.set_ylabel(y1+' - '+y2)
        ax2.plot([-0.5,1.5],[3,12])
        ax2.scatter(pickles[x1]-pickles[x2],pickles[y1]-pickles[y2],c='black',s=5)
        for i in range(len(pickles)):
            ax2.annotate(pickles['name'][i][:-4],xy=(pickles[x1][i]-pickles[x2][i],pickles[y1][i]-pickles[y2][i]),size=10)

        plt.show()
        
