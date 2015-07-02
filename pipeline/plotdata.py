def plotdata(x,y):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.ticker import NullFormatter

	# the random data
	nullfmt   = NullFormatter()         # no labels
	
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02
	
	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	
	# start with a rectangular Figure
	plt.figure(1, figsize=(8,8))
	
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
	xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
	lim = ( int(xymax/binwidth) + 1) * binwidth
	
	#axScatter.set_xlim( (-lim, lim) )
	#axScatter.set_ylim( (-lim, lim) )

	axScatter.set_xlim( (-0.4, 0.6) )
	axScatter.set_ylim( (-0.3, 0.9) )
	
	bins = np.arange(-lim, lim + binwidth, binwidth)
	axHistx.hist(x, bins=bins)
	axHisty.hist(y, bins=bins, orientation='horizontal')
	
	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	
	axScatter.axhline(y=0.,xmin=-2,xmax=2,color='k')
	axScatter.axvline(x=0.,ymin=-2,ymax=2,color='k')

	return plt.show()
