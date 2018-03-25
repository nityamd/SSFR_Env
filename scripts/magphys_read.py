
import numpy as np
import pyfits
import math
from math import log10
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic
from scipy.stats import binned_statistic_2d
from scipy import interpolate
from matplotlib import rc
import astropy as ap
from astropy.cosmology import WMAP7

file = [str(ind2[i]) + '.0' + '.fit' for i in range(len(ind2))]
lines = [open(file[i]).readlines() for i in range(len(ind2))]
print np.shape(lines)
ble = np.vectorize(np.float)
boo = np.asarray([ble(lines[i][10].split()) for i in range(len(ind2))])
#3: TauV 2: mu 4: SSFR 6: Ldust in solar luminosities
mu = boo[:,2]
tauV = boo[:,3]
ssfr = boo[:,4]
log_ssfr = [log10(boo[:,4][i]) for i in range(len(ind2))]
ld = boo[:,6]
mutauv = [mu[i]*tauV[i] for i in range(len(ind2))]
blergh = [(tauV[i] - mu[i]*tauV[i]) for i in range(len(mu))]
