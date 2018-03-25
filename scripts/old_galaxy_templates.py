import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic_2d

data = np.loadtxt("data_Nov_18")

#Just to recap the way that data file is organized is: id, z, c1,c2, opt_flux[7](FNugriz), errors[7], infr_flux[4], errors[4]

c1 = data[:,2]
c2 = data[:,3]

dat = np.asarray([data[i] for i in range(len(c1)) if 0 <= c1[i] <= 8.0 and -2.5 <= c2[i] <= 2.5])
#print np.shape(dat)
index = dat[:,0]
redshift = dat[:,1]
#r-band flux:
r = dat[:,8]

for i in range(len(index)):
    moo = dat[:,8][i]
    for j in range(4,26):
        blah = dat[:,j][i]
        dat[:,j][i] = blah/moo


lol = np.isnan(dat)
lol1 = np.where(lol==True)
lol2 = np.unique(lol1[0])


dats = np.delete(dat,lol2,0)
ind = dats[:,0]
z = dats[:,1]
c1 = dats[:,2]
c2 = dats[:,3]
print len(ind)
goo = binned_statistic_2d(c2,c1,ind,statistic = 'count', bins = (25,25))



b = np.reshape(np.asarray(goo[3]),(len(ind),1))
f = np.append(dats,b, axis = 1)

f = f[f[:,26].argsort()]
a = np.unique(f[:,26])
n = [np.where(f[:,26] == a[i]) for i in range(len(a))]

h = np.array([np.mean(f[n[i]], axis = 0) for i in range(len(a))])
print np.shape(h)
stuff = zip(np.arange(625),np.zeros(625),h[:,4],h[:,11],h[:,5],h[:,12],h[:,6],h[:,13],h[:,7],h[:,14],h[:,8],h[:,15],h[:,9],h[:,16],h[:,10],h[:,17],h[:,18],h[:,22],h[:,19],h[:,23],h[:,20],h[:,24],h[:,21],h[:,25])
color_stuff = zip(np.arange(625),h[:,1],h[:,2],h[:,3])

print np.shape(stuff)

number = np.ravel(goo[0])
ind1 = np.where(number !=0)
number1 = number[ind1]
print 'no of non-zero bins', len(number1)

ind2 = np.where(number1> 10)
number2 = number1[ind2]
print 'no >10', len(number2)

stuff1 = np.asarray(stuff)[ind2]

color_stuff1 = np.asarray(color_stuff)[ind2]
print np.shape(stuff1)

print 'Indices', ind2

print 'Checking',stuff1[:,0]

#let's make a file for magphys
#filename = "mag_input_25_by_25_bins_Nov_18"
#fi = open(filename, 'w')
#for line in stuff1:
#fi.write("  ".join(str(x) for x in line) + "\n")
#fi.close()

#filename1 = "magphys_colors_25_by_25_Nov_18"
#f1 = open(filename1,'w')
#for line in color_stuff1:
#f1.write("  ".join(str(x) for x in line) + "\n")
#f1.close()
