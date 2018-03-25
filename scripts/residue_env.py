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
res = np.asarray([mssfr[i] - uvssfr[i] for i in range(len(mssfr))])

env = [log10(nn[i]) for i in range(len(nn))]

ind = np.where(res<=-5.0)[0]
print b[:,0][ind]

h = binned_statistic(res,nn, statistic = 'mean',bins = 100)
b = [(h[1][i] + h[1][i+1])*0.5 for i in range(100)]
envz = np.log10(h[0])

plt.scatter(res, env)
plt.show()

plt.scatter(b,envz)
plt.show()
