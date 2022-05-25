##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 1d Likelihood profile of $m_A$ with $m_H$ minimized at each point
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc
from iminuit import Minuit

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)
import iminuit
print("iminuit version:", iminuit.__version__)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Values
Scen = 0.15
Ssigma = 0.08
Tcen = 0.27
Tsigma = 0.06
STcorrelation = 0.93

# Output files
outputplot = validation_dir+"mHmA_ST_minuit.pdf"

# Scan ranges
#mA_min = 600
#mA_max = 1400
#mH_min = 600
#mH_max = 1400
mA_min = 100
mA_max = 2100
mH_min = 100
mH_max = 2100

# Minimization ranges
mA_mz_min = 200
mA_mz_max = 2000
mH_mz_min = 200
mH_mz_max = 2000

# Starting values for fit parameters
mH0 = 1000
mA0 = 1000
# Seems like it doesn't matter, even if the values are outside of the limits range

# Fixed value of mHpm
my_mHpm = 1000

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

######################################################################
# Likelihood Calculation
######################################################################

def func(mH, mA):
		z10, z20 = Scen, Tcen
		sig1m, sig1p = Ssigma, Ssigma
		sig2m, sig2p = Tsigma, Tsigma
		p = STcorrelation

		z1, z2 = calc.Scalc(mh = 125, mH = mH, mA = mA, mHpm = my_mHpm, sinba = 1), calc.Tcalc(mh = 125, mH = mH, mA = mA, mHpm = my_mHpm, sinba = 1)

		V1 = sig1p * sig1m
		V1e = sig1p - sig1m
		V2 = sig2p * sig2m
		V2e = sig2p - sig2m
		V1f = V1 + V1e * (z1 - z10)
		V2f = V2 + V2e * (z2 - z20)
		L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
		return L2t

######################################################################
# Minuit scan
######################################################################

print("\n***** performing model fit: migrad, hesse, minos *****")

# Initialize the fit; parameter starting values and limits
m = Minuit(func, mH=mH0, mA=mA0)
m.limits = [(mH_mz_min, mH_mz_max), (mA_mz_min, mA_mz_max)]

# Minimization and error estimation
m.migrad()
m.hesse()   # run covariance estimator
m.minos()

# Display parameter values at the best-fit point
#print("\nbest-fit point:", m.values) 
#print("\nHesse errors:", m.errors)
#print("\nMinos errors:")
#for key, value in list(m.merrors.items()):
#    print(key, value)

my_subtract_min = True

print("\n***** getting the 1-dimensional likelihood profile of mH *****")
xH,yH,rH = m.mnprofile('mH', size=300, bound=(mH_min,mH_max), subtract_min=my_subtract_min)

print("\n***** getting the 1-dimensional likelihood profile of mA *****")
xA,yA,rA = m.mnprofile('mA', size=300, bound=(mA_min,mA_max), subtract_min=my_subtract_min)

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure(figsize=(16,8))

# mHHHHHHHH
ax = fig.add_subplot(121)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)

# Title, labels, color bar...
plt.plot(xH,yH,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([mH_min,mH_max])
plt.xlabel(r'$m_H$[GeV]',fontsize=16)
plt.ylabel(r'$\Delta (-2\log L)$',fontsize=16)
plt.text(200, 7.25, f"1d Likelihood profile of $m_H$ with $m_A$ minimized at each point", fontsize=12)
plt.text(200, 7, f"range of $m_A$ = ({mA_mz_min},{mA_mz_max})", fontsize=12)
plt.text(200, 6.75, fr"$\cos(\beta - \alpha) = 0$, $m_{{H^{{\pm}}}} = {my_mHpm}$ GeV", fontsize=12) # à changer avec my_mHpm
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

# best fit point
#plt.plot([mHmin],[m2logLmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# mAAAAAAAAA
ax = fig.add_subplot(122)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)

# Title, labels, color bar...
plt.plot(xA,yA,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([mA_min,mA_max])
plt.xlabel(r'$m_A$[GeV]',fontsize=16)
plt.ylabel(r'$\Delta (-2\log L)$',fontsize=16)
plt.text(200, 7.25, f"1d Likelihood profile of $m_A$ with $m_H$ minimized at each point", fontsize=12)
plt.text(200, 7, f"range of $m_H$ = ({mH_mz_min},{mH_mz_max})", fontsize=12)
plt.text(200, 6.75, fr"$\cos(\beta - \alpha) = 0$, $m_{{H^{{\pm}}}} = {my_mHpm}$ GeV", fontsize=12) # à changer avec my_mHpm
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

# best fit point
#plt.plot([mAmin],[m2logLmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("result see " + validation_dir)
print("***** done *****\n")
