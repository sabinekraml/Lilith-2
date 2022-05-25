##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 2d likelihood contour on the $m_H$, $m_A$ plane with $m_Hpm$ minimized at each point
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
outputplot = validation_dir+"mHmA_ST_minuit_mHpm.pdf"

# Scan ranges
mA_min = 600
mA_max = 1400
mH_min = 600
mH_max = 1400
mHpm_min = 600
mHpm_max = 1400
#mA_min = 100
#mA_max = 2100
#mH_min = 100
#mH_max = 2100
#mHpm_min = 100
#mHpm_max = 2100

# Minimization ranges # They change everything
mA_mz_min = 46
mA_mz_max = 3000
mH_mz_min = 46
mH_mz_max = 3000
mHpm_mz_min = 999.999
mHpm_mz_max = 1000.0001
#mA_mz_min = 600
#mA_mz_max = 1400
#mH_mz_min = 600
#mH_mz_max = 1400
#mHpm_mz_min = 999.999
#mHpm_mz_max = 1000.0001

# Starting values for fit parameters
#mH0 = 500
#mA0 = 500
#mHpm0 = 500
mH0 = 1017.745
mA0 = 1400
mHpm0 = 1000
#mH0 = 904
#mA0 = 904
#mHpm0 = 1000

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

######################################################################
# Likelihood Calculation
######################################################################

def func(mH, mA, mHpm):
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
#		print("mH, mA, mHpm, L2t = ", mH, mA, mHpm, L2t)
		return L2t

######################################################################
# Minuit scan
######################################################################

print("\n***** performing model fit: migrad, hesse, minos *****")

# Initialize the fit; parameter starting values and limits
m = Minuit(func, mH=mH0, mA=mA0, mHpm=mHpm0)
m.limits = [(mH_mz_min, mH_mz_max), (mA_mz_min, mA_mz_max), (mHpm_mz_min, mHpm_mz_max)]

# Minimization and error estimation
print("start migrad")
m.migrad()
print("\nbest-fit point:", m.values) 
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
#print("start hesse")
#m.hesse()   # run covariance estimator
#print("start minos")
#m.minos()
#print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

# Display parameter values at the best-fit point

#print("\nHesse errors:", m.errors)
#print("\nMinos errors:")
#for key, value in list(m.merrors.items()):
#    print(key, value)

print("\n***** getting the 2-dimensional likelihood profile of mH-mA *****")
#print("test = ", m.mnprofile('mH', size=300, bound=(mH_min,mH_max), subtract_min=True) )
#print("test = ", m.contour('mH','mA', size=100) )
#print("test = ", m.mncontour('mH','mA', size=100) )
ctr_xy = m.mncontour("mH","mA", cl=0.68)
#print(ctr_xy)
#print(ctr_xy[:,0])

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)

# Plot
plt.plot(ctr_xy[:,0],ctr_xy[:,1],color="b",linewidth=2.5)

# best fit point
#plt.plot([mHmin],[mAmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$m_H$[GeV]',fontsize=16)
plt.ylabel(r'$m_A$[GeV]',fontsize=16)
#plt.text(-360, 335, r'$\cos(\beta - \alpha) = 0$, $m_{H^{\pm}} = 1TeV$', fontsize=9)
#plt.text(-360, 285, f"Best points ($\Delta m_H, \Delta m_A$) = ({mHmin:.0f}, {mAmin:.0f})", fontsize=9)
#plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + 1000, mA = mAmin + 1000, mHpm = 1000, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + 1000, mA = mAmin + 1000, mHpm = 1000, sinba = 1):.3f}", fontsize=9)
#plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
#plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "validations/STU/")
print("***** done *****")
