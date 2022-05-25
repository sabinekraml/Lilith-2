##################################################################
#
# The 1- and 2-σ regions allowed by the oblique parameter fit (S, T and U separately) for 2HDM in the alignment limit in the plane of
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
Scen = 0.06
Ssigma = 0.10
Tcen = 0.11
Tsigma = 0.12
Ucen = 0.14
Usigma = 0.09

# Output files
outputS = validation_dir+"mHmA_S_withU.out"
outputT = validation_dir+"mHmA_T_withU.out"
outputU = validation_dir+"mHmA_U_withU.out"
outputplot = validation_dir+"mHmA_S-T-U.pdf"

# Scan ranges
mH_min = -900
mH_max = 900
mA_min = -900
mA_max = 900
#taille = 955
#mH_min = -taille
#mH_max = taille
#mA_min = -taille
#mA_max = taille

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# Likelihood Calculation
######################################################################

def func(X, cen, sig):
    return (X-cen)**2/sig**2

######################################################################
# Scan routine
######################################################################

m2logLminS=10000
m2logLminT=10000
m2logLminU=10000
max=-1

print("***** running scan *****")

fresults = open(outputS, 'w')

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logL = func(calc.Scalc(mh = 125, mH = mH + 1000, mA = mA + 1000, mHpm = 1000, sinba = 1), Scen, Ssigma)
#mH = mH + 1000 and same for mA cause it's delta mH actually so mH-mHpm, with here mHpm = 1000
        if m2logL < m2logLminS:
            m2logLminS = m2logL
            mHminS = mH
            mAminS = mA
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logL + '\n')
fresults.close()

print("minimum at mHS, mAS, -2logL_minS = ", mHminS, mAminS, m2logLminS)

fresults = open(outputT, 'w')

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logL = func(calc.Tcalc(mh = 125, mH = mH + 1000, mA = mA + 1000, mHpm = 1000, sinba = 1), Tcen, Tsigma)
#mH = mH + 1000 and same for mA cause it's delta mH actually so mH-mHpm, with here mHpm = 1000
        if m2logL < m2logLminT:
            m2logLminT = m2logL
            mHminT = mH
            mAminT = mA
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logL + '\n')
fresults.close()

print("minimum at mHT, mAT, -2logL_minT = ", mHminT, mAminT, m2logLminT)

fresults = open(outputU, 'w')

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logL = func(calc.Ucalc(mh = 125, mH = mH + 1000, mA = mA + 1000, mHpm = 1000, sinba = 1), Ucen, Usigma)
#mH = mH + 1000 and same for mA cause it's delta mH actually so mH-mHpm, with here mHpm = 1000
        if m2logL < m2logLminU:
            m2logLminU = m2logL
            mHminU = mH
            mAminU = mA
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logL + '\n')
fresults.close()

print("minimum at mHU, mAU, -2logL_minU = ", mHminU, mAminU, m2logLminU)

print("***** scan finalized *****")

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 6
matplotlib.rcParams['ytick.major.pad'] = 6

fig = plt.figure(figsize=(16,8))


# SSSSSSSSSSSSSSSSSSS
ax = fig.add_subplot(221)

plt.minorticks_on()
plt.tick_params(labelsize=14, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(outputS)

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
plt.plot([mHminS],[mAminS], '+', markersize=15, color = 'black', label = 'best fit')
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\Delta m_H$[GeV]',fontsize=18)
plt.ylabel(r'$\Delta m_A$[GeV]',fontsize=18)
plt.text(-850, 730, r'Regions allowed by S with $\cos(\beta - \alpha) = 0$, $m_{H^{\pm}} = 1TeV$', fontsize=12)
plt.text(-850, 630, f"Best point ($\Delta m_H, \Delta m_A$) = ({mHminS:.0f}, {mAminS:.0f}) with S = {calc.Scalc(mh = 125, mH = mHminS + 1000, mA = mAminS + 1000, mHpm = 1000, sinba = 1):.3f}", fontsize=12)


# TTTTTTTTTTTTTTTTTTTTTTTTTTT
ax = fig.add_subplot(222)

plt.minorticks_on()
plt.tick_params(labelsize=14, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(outputT)

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
plt.plot([mHminT],[mAminT], '+', markersize=15, color = 'black', label = 'best fit')
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\Delta m_H$[GeV]',fontsize=18)
plt.ylabel(r'$\Delta m_A$[GeV]',fontsize=18)
plt.text(-850, 730, r'Regions allowed by T with $\cos(\beta - \alpha) = 0$, $m_{H^{\pm}} = 1TeV$', fontsize=12)
plt.text(-850, 630, f"Best point ($\Delta m_H, \Delta m_A$) = ({mHminT:.0f}, {mAminT:.0f}) with T = {calc.Tcalc(mh = 125, mH = mHminT + 1000, mA = mAminT + 1000, mHpm = 1000, sinba = 1):.3f}", fontsize=12)


# UUUUUUUUUUUUUUUUUUUUUUUU
ax = fig.add_subplot(223)

plt.minorticks_on()
plt.tick_params(labelsize=14, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(outputU)

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
plt.plot([mHminU],[mAminU], '+', markersize=15, color = 'black', label = 'best fit')
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\Delta m_H$[GeV]',fontsize=18)
plt.ylabel(r'$\Delta m_A$[GeV]',fontsize=18)
plt.text(-850, 730, r'Regions allowed by U with $\cos(\beta - \alpha) = 0$, $m_{H^{\pm}} = 1TeV$', fontsize=12)
plt.text(-850, 630, f"Best point ($\Delta m_H, \Delta m_A$) = ({mHminU:.0f}, {mAminU:.0f}) with U = {calc.Ucalc(mh = 125, mH = mHminU + 1000, mA = mAminU + 1000, mHpm = 1000, sinba = 1):.3f}", fontsize=12)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("***** done *****\n")
