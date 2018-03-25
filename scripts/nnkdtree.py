import numpy as np
import pyfits
import sys
import scipy as sp
import scipy.spatial as ss
import math


#sys.argvs used:
#[1] the neighbor upto which we want to go: n

#sys.argvs used:
#[1] aperture size = 0.5 Mpc, 1.0 Mpc, 6.0 Mpc etc..
#[2] velocity projection +- 500 km/s, 1000 km/s etc..

data = np.loadtxt("id_z_colors_radec_coords")
points = zip(data[:,6],data[:,7],data[:,8])
z = data[:,1]
ind = data[:,0]
a = data[:,4]
d = data[:,5]
#KD TREEEEE!!!!
tree  = ss.KDTree(points)

#I want to find the distances and projected distances up to the n'th nearest neighbour
#For starters I'll be doing upto the 3rd neighbour

n = int(float(sys.argv[1]))

neighbours = tree.query(points,k = n+1)

#The distance to the nth neighbour
dist = neighbours[0][:,n]

#The projected distances: Let's first get the index of the n'th neighbour

n_index = neighbours[1][:,n]

dist_proj = []
for i in range(len(z)):
    g =(math.sin(math.radians(d[i])))*(math.sin(math.radians(d[n_index[i]])))
    p =(math.cos(math.radians(d[i])))*(math.cos(math.radians(d[n_index[i]])))*(math.cos(math.radians(a[i] - a[n_index[i]])))
    ang_sep = g+p
    d_transverse = (300000.0/70.0)*(z[i]+z[n_index[i]])*(math.sqrt(abs((1 -ang_sep)/2)))
    d_los =  (300000.0)* abs(z[i] - z[n_index[i]])
    dist_proj.append(d_transverse)


print len(dist), len(dist_proj)

filename = 'dist_' + str(sys.argv[1]) + 'th_neighbour'

z_nn = zip(ind,z,dist,dist_proj)

f = open(filename,'w')
for line in z_nn:
    f.write("  ".join(str(x) for x in line) + "\n")
f.close()
