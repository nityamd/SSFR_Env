import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic_2d
from matplotlib import rc
from matplotlib import colors
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata

#from mpl_toolkits.axes_grid1 import ImageGrid
d = np.loadtxt("big_data_file.txt")
print np.shape(d)
#Big file recap
#Let's make the BIG_FILE
# 0.index 1.opt 2.inf 3.mag_ssfr 4.uv_ssfr 5.mass 6.nn 7.dp3 8.dp5 9.dp10

nsa_sample = np.asarray(d[:,0])

#Big file recap
#Let's make the BIG_FILE
# 0.index 1.opt 2.inf 3.mag_ssfr 4.uv_ssfr 5.mass 6.nn 7.dp3 8.dp5 9.dp10


a = pyfits.open("nsa_v0_1_2.fits")
c = a[1].data
nsa = np.asarray(c["NSAID"])
mag = np.asarray(c["ABSMAG"])
ind = np.in1d(nsa,nsa_sample)
new_nsa = nsa[ind]

mag1 = mag[ind]
rmag = np.asarray(mag[:,4])
rmag1 = np.asarray(mag1[:,4])

print min(rmag), max(rmag)
print min(rmag1), max(rmag1)

lum = np.where(rmag1<=-17.5)[0]
print len(rmag1[lum])

d = d[lum]
print np.shape(d)

youwe = d[:,4]
ind = np.where(youwe >= -16.0)[0]
print len(ind)
b = d[ind]

youwe = b[:,4]
ind = np.where(youwe <= -7.0)[0]
print len(ind)
b = b[ind]
print np.shape(b)




opt = b[:,1]
inf = b[:,2]
mssfr = b[:,3]
uvssfr = b[:,4]
glug = np.vectorize(log10)
mass = glug(b[:,5])
nn = b[:,6]
dp3 = b[:,7]
dp5 = b[:,8]
dp10 = b[:,9]

#ro3 = [log10(abs(3.0/(3.14159*((float(dp3[i]))**2.0)))) for i in range(len(dp3))]
#ro5 = [log10(abs(5.0/(3.14159*((float(dp5[i]))**2.0)))) for i in range(len(dp5))]
#ro10 = [log10(abs(10.0/(3.14159*((float(dp10[i]))**2.0)))) for i in range(len(dp10))]

def oden(a):
    b = np.log10((70.0*(float(a) - 1.0))/(3.14159*(0.5**2.0)*2000.0))
    return b
boo = np.vectorize(oden)


levels = np.linspace(-14.0,-8.0,11)

h = binned_statistic_2d(inf,opt,inf,statistic = 'count',bins = (25,25))
num = binned_statistic_2d(inf,opt,mssfr,statistic = 'count',bins = (25,25))
env = binned_statistic_2d(inf,opt,nn,statistic = 'mean',bins = (25,25))
ms = binned_statistic_2d(inf,opt,mssfr,statistic = 'mean',bins = (25,25))
uvs = binned_statistic_2d(inf,opt,uvssfr,statistic = 'mean',bins = (25,25))

number = np.ma.array(h[0], mask = h[0] == 0)
number = np.log10(number)

yedges = h[1]
xedges = h[2]
xbins = [0.5*(h[1][i] + h[1][i+1]) for i in range(25)]
ybins = [0.5*(h[2][i] + h[2][i+1]) for i in range(25)]
extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]

plt.pcolormesh(yedges,xedges,(np.transpose(number)),cmap = plt.cm.binary)
plt.colorbar()
plt.show()


def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    asp = (abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    return asp

f = plt.figure(figsize=(8,8)) # Notice the equal aspect ratio
ax = [f.add_subplot(2,2,i+1) for i in range(4)]
f.subplots_adjust(wspace=0, hspace = 0)



ax1 = ax[0]
ax2 = ax[1]
ax3 = ax[2]
ax4 = ax[3]



cs1 = ax1.pcolormesh(yedges,xedges,(np.transpose(number)),cmap = plt.cm.binary)
ax1.set_xlim(-2.5,2.5)
ax1.set_ylim(0,8.0)
ax1.xaxis.set_ticklabels([])

#cmap1 = plt.cm.binary
##bounds=[0,0.5,1.0,1.5,2.0,2.5,3.0]
##norm = plt.cm.BoundaryNorm(bounds, cmap1.N)
#cs1 = ax1.imshow(np.flipud(np.transpose(number)),interpolation = 'nearest',extent = [-2.5,2.5,0,8.0],cmap=plt.cm.binary)
#ax1.set_xlim(-2.5,2.5)
#ax1.set_ylim(0,8.0)
##ax1.set_title('Number Distribution', fontsize = 12)
#ax1.xaxis.set_ticklabels([])
#asp = forceAspect(ax1,aspect=1)
#ax1.set_aspect(asp)
#ax1.set_yticks([0,1.0,2.0,3.0,4.0,5.0,6.0,7.0])
##wow = np.ma.array(h[0],mask = h[0]==0)

cs2 = ax2.contourf(xbins, ybins, np.transpose(boo(env[0])),10,corner_mask = False,extent = [-2.5,2.5,0,8.0], cmap=plt.cm.bone)
cs2b = ax2.contour(cs2,colors = 'k')
ax2.clabel(cs2b, colors = 'k', fontsize = 8)
ax2.set_xlim(-2.5,2.5)
ax2.set_ylim(0,8.0)
ax2.yaxis.set_ticklabels([])
ax2.xaxis.set_ticklabels([])
#ax2.set_aspect(asp)
#forceAspect(ax2,aspect=1)
#plt.title("MAGPHYS SSFR")

cs3 = ax3.contourf(xbins,ybins, np.transpose(ms[0]),levels = levels,corner_mask = False,extent = [-2.5,2.5,0,8.0],cmap=plt.cm.Blues)
cs3b = ax3.contour(cs3,colors = 'k')
ax3.clabel(cs3b, colors = 'k', fontsize = 8)
ax3.set_xlim(-2.5,2.5)
ax3.set_ylim(0,8.0)
#ax3.set_aspect(asp)
ax3.set_yticks([0,1.0,2.0,3.0,4.0,5.0,6.0,7.0])
#forceAspect(ax3,aspect=1)

cs4 = ax4.contourf(xbins, ybins, np.transpose(uvs[0]), levels = levels,corner_mask = False, extent = [-2.5,2.5,0,8.0],cmap=plt.cm.Blues)
cs4b = ax4.contour(cs4,colors = 'k')
ax4.clabel(cs4b, colors = 'k', fontsize = 8)
ax4.set_xlim(-2.5,2.5)
ax4.set_ylim(0,8.0)
ax4.yaxis.set_ticklabels([])
#ax4.set_aspect(asp)
#forceAspect(ax4,aspect=1)


f.text(0.495, 0.06, '[W1] - [W3]', ha='center', va='center', fontsize = 18)
f.text(0.06, 0.5, 'N-r', ha='center', va='center', rotation='vertical', fontsize = 18)


cbar_ax1 = f.add_axes([0.14, 0.91, 0.33, 0.015])
cb1 = f.colorbar(cs1,orientation = 'horizontal', cax=cbar_ax1)
cb1.ax.tick_params(labelsize=7)
#cb1.ax.xaxis.set_ticks([0,0.5,1.0,1.5,2.0,2.5,3.0])

cb1.ax.xaxis.set_ticks_position('top')
cbar_ax2 = f.add_axes([0.56, 0.91, 0.33, 0.015])
cb2 = f.colorbar(cs2,orientation = 'horizontal', cax=cbar_ax2)
cb2.ax.xaxis.set_ticks_position('top')
cb2.ax.tick_params(labelsize=7)
cbar_ax3 = f.add_axes([0.14, 0.02, 0.33, 0.015])
cb3 = f.colorbar(cs3,orientation = 'horizontal', cax=cbar_ax3)
cb3.ax.tick_params(labelsize=7)
cbar_ax4 = f.add_axes([0.56, 0.02, 0.33, 0.015])
cb4 = f.colorbar(cs4,orientation = 'horizontal', cax=cbar_ax4)
cb4.ax.tick_params(labelsize=7)
#plt.tight_layout()
#f.subplots_adjust(wspace=0, hspace=0)
plt.show()
