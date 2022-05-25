##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 1d Likelihood profile of each mass with the 2 other minimized at each point
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
outputplot = validation_dir+"mHmAmHpm_ST_minuit.pdf"

# Scan ranges
#mA_min = 600
#mA_max = 1400
#mH_min = 600
#mH_max = 1400
#mHpm_min = 600
#mHpm_max = 1400
mA_min = 200
mA_max = 2000
mH_min = 200
mH_max = 2000
mHpm_min = 200
mHpm_max = 2000

# Minimization ranges # They change everything
mA_mz_min = 200
mA_mz_max = 2000
mH_mz_min = 200
mH_mz_max = 2000
mHpm_mz_min = 200
mHpm_mz_max = 2000

# Starting values for fit parameters
mH0 = 500
mA0 = 500
mHpm0 = 500

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

######################################################################
# Likelihood Calculation
######################################################################

def func(mH, mA, mHpm):
#	 print("mH, mA, mHpm = ", mH, mA, mHpm)
		z10, z20 = Scen, Tcen
		sig1m, sig1p = Ssigma, Ssigma
		sig2m, sig2p = Tsigma, Tsigma
		p = STcorrelation

		z1, z2 = calc.Scalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1), calc.Tcalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1)

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
m = Minuit(func, mH=mH0, mA=mA0, mHpm=mHpm0)
m.limits = [(mH_mz_min, mH_mz_max), (mA_mz_min, mA_mz_max), (mHpm_mz_min, mHpm_mz_max)]

# Minimization and error estimation
#m.migrad()
#m.hesse()   # run covariance estimator
#m.minos()


# Display parameter values at the best-fit point
#print("\nbest-fit point:", m.values) 
#print("\nHesse errors:", m.errors)
#print("\nMinos errors:")
#for key, value in list(m.merrors.items()):
#    print(key, value)

my_subtract_min = True

print("\n***** getting the 1-dimensional likelihood profile of mH *****")
# Profiling over mH
xH,yH,rH = m.mnprofile('mH', size=300, bound=(mH_min,mH_max), subtract_min=my_subtract_min)
#print("rH = ", rH)

print("\n***** getting the 1-dimensional likelihood profile of mA *****")
# Profiling over mA
xA,yA,rA = m.mnprofile('mA', size=300, bound=(mA_min,mA_max), subtract_min=my_subtract_min)
print("\n***** getting the 1-dimensional likelihood profile of mHpm *****")
# Profiling over mHpm
xHpm,yHpm,rHpm = m.mnprofile('mHpm', size=300, bound=(mHpm_min,mHpm_max), subtract_min=my_subtract_min)

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure(figsize=(16,8))

# mHHHHHHHH
ax = fig.add_subplot(221)

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
plt.text(250, 7.1, r'1d Likelihood profile of $m_H$ with $\cos(\beta - \alpha) = 0$', fontsize=13)
plt.text(250, 6.6, r'$m_A$ and $m_{H^{\pm}}$ minimized at each point', fontsize=13)
plt.text(250, 6.1, f"range of $m_A$, $m_{{H^{{\pm}}}}$ = ({mA_mz_min},{mA_mz_max}),({mHpm_mz_min},{mHpm_mz_max})", fontsize=13) # to update with $m_{H^{\pm}}$
plt.axhline(y=2.3,color='k',ls='dashed') # Check the values
plt.axhline(y=5.99,color='k',ls='dashed')

# best fit point
#plt.plot([mHmin],[m2logLmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# mAAAAAAAAA
ax = fig.add_subplot(222)

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
plt.text(250, 7.1, r'1d Likelihood profile of $m_A$ with $\cos(\beta - \alpha) = 0$', fontsize=13)
plt.text(250, 6.6, r'$m_H$ and $m_{H^{\pm}}$ minimized at each point', fontsize=13)
plt.text(250, 6.1, fr"range of $m_H$, $m_{{H^{{\pm}}}}$ = ({mH_mz_min},{mH_mz_max}),({mHpm_mz_min},{mHpm_mz_max})", fontsize=13) # to update with $m_{H^{\pm}}$
plt.axhline(y=2.3,color='k',ls='dashed')
plt.axhline(y=5.99,color='k',ls='dashed')

# best fit point
#plt.plot([mAmin],[m2logLmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# mHHHHHHHHHHHpm
ax = fig.add_subplot(223)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)

# Title, labels, color bar...
plt.plot(xHpm,yHpm,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([mHpm_min,mHpm_max])
plt.xlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=16)
plt.ylabel(r'$\Delta (-2\log L)$',fontsize=16)
plt.text(250, 7.1, r'1d Likelihood profile of $m_{H^{\pm}}$ with $\cos(\beta - \alpha) = 0$', fontsize=13)
plt.text(250, 6.6, r'$m_A$ and $m_H$ minimized at each point', fontsize=13)
plt.text(250, 6.1, f"range of $m_A$, $m_H$ = ({mA_mz_min},{mA_mz_max}),({mH_mz_min},{mH_mz_max})", fontsize=13)
plt.axhline(y=2.3,color='k',ls='dashed')
plt.axhline(y=5.99,color='k',ls='dashed')

# best fit point
#plt.plot([mHpmmin],[m2logLmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("result see " + validation_dir)
print("***** done *****\n")
