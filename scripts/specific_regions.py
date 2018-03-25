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
b = np.loadtxt("big_data_file")
#Big file recap
#Let's make the BIG_FILE
# 0.index 1.opt 2.inf 3.mag_ssfr 4.uv_ssfr 5.mass 6.nn 7.dp3 8.dp5 9.dp10
nsaid = b[:,0]
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
#
#def oden(a):
#    b = np.log10((70.0*(float(a) - 1.0))/(3.14159*(0.5**2.0)*2000.0))
#    return b
#boo = np.vectorize(oden)
#levels = np.linspace(-14.0,-8.0,11)
res = [mssfr[i] - uvssfr[i] for i in range(len(mssfr))]

h = binned_statistic_2d(inf,opt,inf,statistic = 'count',bins = (25,25))
#num = binned_statistic_2d(inf,opt,mssfr,statistic = 'count',bins = (25,25))
#env = binned_statistic_2d(inf,opt,nn,statistic = 'mean',bins = (25,25))
#ms = binned_statistic_2d(inf,opt,mssfr,statistic = 'mean',bins = (25,25))
#uvs = binned_statistic_2d(inf,opt,uvssfr,statistic = 'mean',bins = (25,25))
#residue = binned_statistic_2d(inf,opt,res,statistic = 'mean',bins = (25,25))

yedges = h[1]
xedges = h[2]

k = np.unique(h[3])
#print k

#for i in range(len(k)):
#    ind = np.where(h[3]==k[i])[0]
#    ble = inf[ind]
#    mle = opt[ind]
#    print k[i], min(ble), max(ble), min(mle), max(mle)


indices = [400,401,402,427,428,429,430,454,455,481,482,508]

ids =[]

for i in range(len(indices)):
    ind = np.asarray(np.where(h[3]==indices[i])[0])
    ids.append(ind)

lst  = np.concatenate((ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7],ids[8],ids[9],ids[10],ids[11]))

print lst

bnew = b[lst]

print bnew

filename = 'funky'

f = open(filename,'w')

for line in bnew:
    f.write("  ".join(str(x) for x in line) + "\n")

f.close()
