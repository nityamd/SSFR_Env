import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic_2d
from matplotlib import rc
import pickle

#Steps: Starting: (145155)
#Make a volume limited sample; (-18.5> mr > -24.5); (95738, 95638)
#Remove weird Galex fluxes;(77547,)(76769,)
#Remove weird Optical fluxes; (76765,)(76763,)(76752,)(76737,)(76710)
#Remove weird WISE fluxes; 76693, 76683, 76645, 76607
#(N-r) cut: 7.5>N-r; (75536)   N-r>=0; (75533)
#[W1]-[W3] cut: (-3.0,3.0); (75505, 75476)


fi = pyfits.open("nsa_wise_v0_2_1.fits")
dat = fi[1].data
d = np.asarray(dat)
#-----------------------------------------------------------------------------------
#defining the volume limited sample
#r-band magnitudes
mr = d['ABSMAG'][:,4]
indr1 = np.where(mr<=-18.5)[0]
d = d[indr1]
mr = d['ABSMAG'][:,4]
indr2 = np.where(mr>=-24.5)[0]
d = d[indr2]



#make an environment defining population file
#Things in the file: 0. NSAID; 1. z; 2. RA; 3. Dec;
stuff = np.column_stack((d['NSAID'],d['Z'],d['RA'],d['DEC']))
#filename = 'edp'
#f = open(filename,'w')
#for line in stuff:
#    f.write("  ".join(str(x) for x in line) + "\n")
#f.close()

#Let's also make a stuctured array of all the galaxies in EDP with their other properties too;
# i.e. Pickling time
afile0 = open(r'EDP.pkl', 'wb')
pickle.dump(d, afile0)
afile0.close()


#------------------------------------------------------------------------------------
#removing bad galex stuff; then optical;
flux = d['NMGY']
fuv = flux[:,0]
ind1 = np.where(fuv>0)[0]
d = d[ind1]
ind2 = np.where(d['NMGY'][:,1]>0)[0]
d = d[ind2]
ind3 = np.where(d['NMGY'][:,2]>0)[0]
d = d[ind3]
ind4 = np.where(d['NMGY'][:,3]>0)[0]
d = d[ind4]
ind5 = np.where(d['NMGY'][:,4]>0)[0]
d = d[ind5]
ind6 = np.where(d['NMGY'][:,5]>0)[0]
d = d[ind6]
ind7 = np.where(d['NMGY'][:,6]>0)[0]
d = d[ind7]

#WISE fluxes; Replacing zero fluxes with forced and then filtering out the non-zero fluxes

#W1
w1 = np.asarray(d['W1_NMGY'])
w1f = np.asarray(d['W1_FORCED_NMGY'])
indw1 = np.where(w1==0)[0]
lug = np.arange(len(w1))
indz = lug[indw1]
np.put(w1,indz,w1f[indw1])
indw11 = np.where(w1>0)[0]
d = d[indw11]
w1 = w1[indw11]
d['W1_NMGY'] = w1
#W2
w2 = np.asarray(d['W2_NMGY'])
w2f = np.asarray(d['W2_FORCED_NMGY'])
indw2 = np.where(w2==0)[0]
lug2 = np.arange(len(w2))
indz2 = lug2[indw2]
np.put(w2,indz2,w2f[indw2])
indw22 = np.where(w2>0)[0]
d = d[indw22]
w2 = w2[indw22]
d['W2_NMGY'] = w2
#W3
w3 = np.asarray(d['W3_NMGY'])
w3f = np.asarray(d['W3_FORCED_NMGY'])
indw3 = np.where(w3==0)[0]
lug3 = np.arange(len(w3))
indz3 = lug3[indw3]
np.put(w3,indz3,w3f[indw3])
indw33 = np.where(w3>0)[0]
d = d[indw33]
w3 = w3[indw33]
d['W3_NMGY'] = w3
#W4
w4 = np.asarray(d['W4_NMGY'])
w4f = np.asarray(d['W4_FORCED_NMGY'])
indw4 = np.where(w4==0)[0]
lug4 = np.arange(len(w4))
indz4 = lug4[indw4]
np.put(w4,indz4,w4f[indw4])
indw44 = np.where(w4>0)[0]
d = d[indw44]
w4 = w4[indw44]
d['W4_NMGY'] = w4


#------------------------------------------------------------------------------------
#Restricting Optical colors to 0<= N-r <= 7.5
opt = np.asarray([d['ABSMAG'][:,1][i] - d['ABSMAG'][:,4][i] for i in range(len(d['NSAID']))])
inf = np.asarray([d['W1MAG'][i] - d['W3MAG'][i] for i in range(len(d['NSAID']))])

indo1 = np.where(opt<=7.5)[0]
do = d[indo1]
opto = opt[indo1]
indo2 = np.where(opto>=0)[0]
do = do[indo2]

#Restricting IR colors to -3.0 <= [W1]-[W3] <= 3.0
inf = np.asarray([do['W1MAG'][i] - do['W3MAG'][i] for i in range(len(do['NSAID']))])
indf1 = np.where(inf>=-3.0)[0]
infi = inf[indf1]
di = do[indf1]
indf2 = np.where(infi<=3.0)[0]
di=di[indf2]

#------------------------------------------------------------------------------------
#PICKLING THE ND ARRAY
print(np.shape(di))

#opt = np.asarray([di['ABSMAG'][:,1][i] - di['ABSMAG'][:,4][i] for i in range(len(d['NSAID']))])
#inf = np.asarray([di['W1MAG'][i] - do['W3MAG'][i] for i in range(len(do['NSAID']))])

#afile = open(r'sfr_population.pkl', 'wb')
#pickle.dump(di, afile)
#afile.close()

