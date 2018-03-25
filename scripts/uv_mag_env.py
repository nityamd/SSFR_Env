import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic
from scipy.stats import binned_statistic_2d
from matplotlib import rc
from matplotlib import colors
from matplotlib.colors import LogNorm

d = np.loadtxt("big_data_file.txt")
print np.shape(d)
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
#oden = np.array([ log10((70.0*float(nn[i]))/(3.14159*(0.5**2.0)*2000.0)) for i in range(len(nn))])

#def oden(a):
#    b = np.log10((70.0*(float(a) - 1.0))/(3.14159*(0.5**2.0)*2000.0))
#    return b
#boo = np.vectorize(oden)
#



h = binned_statistic_2d(mssfr,uvssfr,nn,statistic = 'mean',bins = (10,10))
yedges = h[1]
xedges = h[2]


xbins = [0.5*(h[1][i] + h[1][i+1]) for i in range(10)]
ybins = [0.5*(h[2][i] + h[2][i+1]) for i in range(10)]


g = binned_statistic_2d(mssfr,uvssfr,nn,statistic = 'count',bins = (10,10))
yedges = g[1]
xedges = g[2]

#print g[0]

xbins = [0.5*(g[1][i] + g[1][i+1]) for i in range(10)]
ybins = [0.5*(g[2][i] + g[2][i+1]) for i in range(10)]

#extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]


number = np.ma.array(h[0], mask = h[0] <= 10)

gee = g[0]
ech = h[0]
#gee = np.ma.array(g[0], mask = g[0] <= 50)
#ech = np.ma.array(h[0], mask = g[0] <= 50)

x = [-7.0,-8.5,-13.5, -16.5]
y = [-7.0,-8.5,-13.5, -16.5]


#plt.imshow(np.flipud(np.transpose(np.log10(gee))), aspect = 1, extent = extent, cmap = plt.cm.binary_r)
#plt.plot(x,y, color = 'k')
#plt.xlim(-14.5,-7.5)
#plt.ylim(-15.0,-7.0)
#plt.title("Distribution of Field Galaxies (Log10(number) color-coded)")
#plt.xlabel("MAGPHYS SSFR")
#plt.ylabel("UV SSFR")
#plt.colorbar()
#plt.show()

f, (ax1, ax2) = plt.subplots(1, 2, sharex = True, sharey=True)

bounds = np.linspace(0, 1.6, 21)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

cs1 = ax1.pcolormesh(yedges,xedges, np.ma.masked_invalid(np.transpose(np.log10(ech))),norm = norm,cmap=plt.cm.jet)
ax1.plot(x,y)
ax1.set_xlim(-15.5,-7.5)
ax1.set_ylim(-16.5,-7.0)
ax1.set_title("log10(No. of Nearest Neighbors)")
cs2 = ax2.pcolormesh(yedges,xedges, np.transpose(np.ma.masked_invalid(np.log10(gee))),cmap = plt.cm.binary_r)
ax2.plot(x,y)
ax2.set_xlim(-15.5,-7.5)
ax2.set_ylim(-16.5,-7.0)
ax2.set_title("log10(Number of galaxies)")
cbar_ax1 = f.add_axes([0.14, 0.91, 0.33, 0.015])
cb1 = f.colorbar(cs1,orientation = 'horizontal', cax=cbar_ax1)
cb1.ax.tick_params(labelsize=7)
cb1.ax.xaxis.set_ticks_position('top')
cbar_ax2 = f.add_axes([0.56, 0.91, 0.33, 0.015])
cb2 = f.colorbar(cs2,orientation = 'horizontal', cax=cbar_ax2)
cb2.ax.xaxis.set_ticks_position('top')
cb2.ax.tick_params(labelsize=7)
f.text(0.465, 0.04, 'MAGPHYS SSFR', ha='center', va='center', fontsize = 17)
f.text(0.06, 0.5, 'UV SSFR', ha='center', va='center', rotation='vertical', fontsize = 17)
plt.show()
