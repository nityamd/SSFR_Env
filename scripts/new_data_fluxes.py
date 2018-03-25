import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic_2d


#This is for the updated catalogue v0_2_1

fi = pyfits.open("nsa_wise_v0_2_1.fits")
dat = fi[1].data


#Let's just first remove all the galaxies with "zeroes" for any of the flux values
#Note that we're just using the k-corrected fluxes for NSA data and apparent magnitudes for WISE

w1 = dat["W1_NMGY"]
w1f = dat["W1_FORCED_NMGY"]
ind1 = np.where(w1==0)[0]
print 'w1', len(ind1)
for i in range(len(ind1)):
    w1[ind1[i]] = w1f[ind1[i]]
print 'new_w1', len(np.where(w1==0)[0])

w2 = dat["W2_NMGY"]
w2f = dat["W2_FORCED_NMGY"]
ind2 = np.where(w2==0)[0]
print 'w2', len(ind2)
for i in range(len(ind2)):
    w2[ind2[i]] = w2f[ind2[i]]
print 'new_w2', len(np.where(w2==0)[0])

w3 = dat["W3_NMGY"]
w3f = dat["W3_FORCED_NMGY"]
ind3 = np.where(w3==0)[0]
print 'w3', len(ind3)
for i in range(len(ind3)):
    w3[ind3[i]] = w3f[ind3[i]]
print 'new_w3', len(np.where(w3==0)[0])

w4 = dat["W4_NMGY"]
w4f = dat["W4_FORCED_NMGY"]
ind4 = np.where(w4==0)[0]
print 'w4', len(ind4)
for i in range(len(ind4)):
    w4[ind4[i]] = w4f[ind4[i]]
print 'new_w4', len(np.where(w4==0)[0])

combined_data = np.column_stack((dat["NMGY"],w1,w2,w3,w4,dat["NMGY_IVAR"]))

#the indexes that remain after removing the "zeros" shall be called inds
inds = np.all(combined_data != 0, axis = 1)

print len(inds)

#NSA OPTICAL STUFF
#Let's first get the k-corrected fluxes in janskys
flux = dat["NMGY"][inds]
kc = dat["KCORRECT"][inds]

def jansky(a,b):
    ma = a*3631*(10.0**(-9.0))*(10**(b/(-2.5)))
    return ma
jansky = np.vectorize(jansky)
flux = jansky(flux,kc)
print len(flux)
print np.shape(flux)
n = len(flux)

#Flux inverse variance all in units of nanomaggies^-2
flerr = dat["NMGY_IVAR"][inds]
def jansky_err(a,b):
    ma_err = (a**(-0.5))*3631*(10.0**(-9.0))*(10**(b/(-2.5)))
    return ma_err
jansky_err = np.vectorize(jansky_err)
flerr = jansky_err(flerr,kc)

#Let's get the colors,NSA Ids, redshifts, etc
amag = dat["ABSMAG"][inds]
z = dat["Z"][inds]
id = dat["NSAID"][inds]
#N-r
c1 = [(amag[:,1][i] - amag[:,4][i]) for i in range(n)]

#WISE INFRARED STUFF
#WISE fluxes in nanomaggies
wf = np.column_stack((w1[inds],w2[inds],w3[inds],w4[inds]))
#WISE fluxes errors in nanomaggies
wf_err = np.column_stack((dat["W1SIGM_NMGY"][inds],dat["W2SIGM_NMGY"][inds],dat["W3SIGM_NMGY"][inds],dat["W4SIGM_NMGY"][inds]))
#WISE magnitudes
wmags = np.column_stack((dat["W1MAG"][inds],dat["W2MAG"][inds],dat["W3MAG"][inds],dat["W4MAG"][inds]))
#W1-W3
c2 = [(wmags[:,0][i] - wmags[:,2][i]) for i in range(n)]

def wjansky(a):
    ma = 3631*(10**(-9.0))*a
    return ma
wjansky =  np.vectorize(wjansky)
#WISE fluxes and flux_errors in Janksy:
wflux = wjansky(wf)
wflux_err = wjansky(wf_err)


stuff = np.column_stack((id,z,c1,c2,flux,flerr,wflux,wflux_err))
#print np.shape(stuff)
#filename2 = "data_Nov_18"
#f2 = open(filename2,'w')
#for line in stuff:
#f2.write("  ".join(str(x) for x in line) + "\n")
#f2.close()



#For comparisons with k-corrected stellar masses and getting UV fluxes so we can do the Salim et. al. comparison as well
mass = dat["MASS"][inds]
ha = dat["HAFLUX"][inds]
ha_err = dat["HAFLUXERR"][inds]

#Order: Index, redshift, N-r, W1-W3,mass,FUV,FUV_err,NUV,NUV_err,F_absmag,N_absmag
stuff1 = np.column_stack((id,z,c1,c2,mass,flux[:,0],flerr[:,0],flux[:,1],flerr[:,1],amag[:,0],amag[:,1]))

c1 = stuff1[:,2]
c2 = stuff1[:,3]
dato = np.asarray([stuff1[i] for i in range(len(c1)) if 0 <= c1[i] <= 8.0 and -2.5 <= c2[i] <= 2.5])


print np.shape(dato)
filename3 = "mass_colors_uv_fluxes"
f3 = open(filename3,'w')
for line in dato:
    f3.write("  ".join(str(x) for x in line) + "\n")
f3.close()
