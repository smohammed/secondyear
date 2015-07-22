import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#=======================================================================
# Show value of z axis when plotting
#=======================================================================

"""
Show how to modify the coordinate formatter to report the image "z"
value of the nearest pixel given x and y
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(im5bg,vmin=0,vmax=0.1,origin='lower',interpolation='nearest',aspect='auto',cmap=cm.gray)

numrows, numcols = im5bg.shape
def format_coord(x, y):
    col = int(x+0.5)
    row = int(y+0.5)
    if col>=0 and col<numcols and row>=0 and row<numrows:
        z = im5bg[row,col]
        return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
    else:
        return 'x=%1.4f, y=%1.4f'%(x, y)

ax.format_coord = format_coord
plt.show()
