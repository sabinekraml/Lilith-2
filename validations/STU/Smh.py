##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STUc as STUc

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
Scen = 0.15
Ssigma = 0.08
Tcen = 0.27
Tsigma = 0.06
STcorrelation = 0.93


# Output files
output = validation_dir+"Smh.out"
outputplot = validation_dir+"Smh.pdf"

# Scan ranges
T_min = -0.5
T_max = 0.5
mHpm_min = 0
mHpm_max = 200

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

#def func(X, cen1, cen2, sig1p, sig1m, sig2p, sig2m, p):
#    z1, z2 = X[0], X[1]
#    z10, z20 = cen1, cen2

#    V1 = sig1p * sig1m
#    V1e = sig1p - sig1m
#    V2 = sig2p * sig2m
#    V2e = sig2p - sig2m
#    V1f = V1 + V1e * (z1 - z10)
#    V2f = V2 + V2e * (z2 - z20)
#    L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#    return L2t

######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for mHpm in np.linspace(mHpm_min, mHpm_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen, Tcen, Ssigma, Ssigma, Tsigma, Tsigma, STcorrelation)
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            Smin = S
            Tmin = T
        fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at S, T, -2logL_min = ", Smin, Tmin, m2logLmin)


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

ax.set_aspect((S_max-S_min)/(T_max-T_min))

# best fit point
plt.plot([Smin],[Tmin], '+', markersize=8, color = 'white', label = 'CDF best fit')

# Standard Model 
plt.plot([0], [0], '*',markersize=8, color = 'black', label="SM prediction")
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'S',fontsize=16)
plt.ylabel(r'$T$',fontsize=16)
plt.text(-0.15, 0.45, r'fixed U = 0', fontsize=11)

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/validations/STU/")
print("***** done *****")

