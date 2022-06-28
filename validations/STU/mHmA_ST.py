##################################################################
#
# The 1- and 2-σ regions allowed by the oblique parameter fit (S and T combined) for 2HDM in the alignment limit in the plane of
# ∆mA versus ∆mH, where ∆mA = mA − mH± and ∆mH = mH − mH± 
# using best fit from Table 3 : https://arxiv.org/pdf/2204.03796.pdf
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Values
#Scen = 0.15
#Ssigma = 0.08
#Tcen = 0.27
#Tsigma = 0.06
#STcorrelation = 0.93

Scen = 0.05
Ssigma = 0.08
Tcen = 0.09
Tsigma = 0.07
STcorrelation = 0.92

# Output files
#output = validation_dir+"mHmA_ST_CDF.out"
#outputplot = validation_dir+"mHmA_ST_CDF.pdf"
output = validation_dir+"mHmA_ST_PDG.out"
outputplot = validation_dir+"mHmA_ST_PDG.pdf"

# Scan ranges
mH_min = -400
mH_max = 400
mA_min = -400
mA_max = 400

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

# Prepare output
fresults = open(output, 'w')

######################################################################
# Likelihood Calculation
######################################################################

def func(X, cen1, cen2, sig1p, sig1m, sig2p, sig2m, p):
    z1, z2 = X[0], X[1]
    z10, z20 = cen1, cen2

    V1 = sig1p * sig1m
    V1e = sig1p - sig1m
    V2 = sig2p * sig2m
    V2e = sig2p - sig2m
    V1f = V1 + V1e * (z1 - z10)
    V2f = V2 + V2e * (z2 - z20)
    L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
    return L2t

######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

my_mHpm = 1000

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logL = func([calc.Scalc(mh = 125, mH = mH + my_mHpm, mA = mA + my_mHpm, mHpm = my_mHpm, sinba = 1) , calc.Tcalc(mh = 125, mH = mH + my_mHpm, mA = mA + my_mHpm, mHpm = my_mHpm, sinba = 1) ], Scen, Tcen, Ssigma, Ssigma, Tsigma, Tsigma, STcorrelation)
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            mHmin = mH
            mAmin = mA
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at mH, mA, -2logL_min = ", mHmin, mAmin, m2logLmin)


######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
matplotlib.rcParams['axes.linewidth'] = 2

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.yticks(np.linspace(-400, 400, 5))

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())


# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")


# Plotting the 68% and 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

# best fit point
plt.plot([mHmin],[mAmin], '+', markersize=8, color = 'black', label = 'best fit')
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\Delta m_H$[GeV]',fontsize=16)
plt.ylabel(r'$\Delta m_A$[GeV]',fontsize=16)
plt.text(-360, 335, r'$\cos(\beta - \alpha) = 0$, $m_{H^{\pm}} = 1TeV$', fontsize=9)
plt.text(-360, 285, f"Best points ($\Delta m_H, \Delta m_A$) = ({mHmin:.0f}, {mAmin:.0f})", fontsize=9)
plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}", fontsize=9)
#plt.text(-360, 205, f"CDF (2022)", fontsize=9)
plt.text(-360, 205, f"PDG (2021)", fontsize=9)
#plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
#plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "validations/STU/")
print("***** done *****")
