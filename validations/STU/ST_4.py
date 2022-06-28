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

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))+"/"
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
Scen_PDG = 0.05
Ssigma_PDG = 0.08
Tcen_PDG = 0.09
Tsigma_PDG = 0.07
STcorrelation_PDG = 0.92

Scen_CDF = 0.15
Ssigma_CDF = 0.08
Tcen_CDF = 0.27
Tsigma_CDF = 0.06
STcorrelation_CDF = 0.93

Scen_CDF_std = 0.100
Ssigma_CDF_std = 0.073
Tcen_CDF_std = 0.202
Tsigma_CDF_std = 0.056
STcorrelation_CDF_std = 0.93

Scen_CDF_cons = 0.086
Ssigma_CDF_cons = 0.077
Tcen_CDF_cons = 0.177
Tsigma_CDF_cons = 0.070
STcorrelation_CDF_cons = 0.89

# Output files
output = validation_dir+"ST_4.out"
outputplot = validation_dir+"ST_2.pdf"

# Scan ranges
S_min = -0.2
S_max = 0.4
T_min = -0.1
T_max = 0.5

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

for S in np.linspace(S_min, S_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen_PDG, Tcen_PDG, Ssigma_PDG, Ssigma_PDG, Tsigma_PDG, Tsigma_PDG, STcorrelation_PDG)
        fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f     '%m2logL + '\n')

print("PDG done")

for S in np.linspace(S_min, S_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen_CDF, Tcen_CDF, Ssigma_CDF, Ssigma_CDF, Tsigma_CDF, Tsigma_CDF, STcorrelation_CDF)
        fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f     '%m2logL + '\n')

print("CDF done")

for S in np.linspace(S_min, S_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen_CDF_std, Tcen_CDF_std, Ssigma_CDF_std, Ssigma_CDF_std, Tsigma_CDF_std, Tsigma_CDF_std, STcorrelation_CDF_std)
        fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f     '%m2logL + '\n')

print("CDF standard done")

for S in np.linspace(S_min, S_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen_CDF_cons, Tcen_CDF_cons, Ssigma_CDF_cons, Ssigma_CDF_cons, Tsigma_CDF_cons, Tsigma_CDF_cons, STcorrelation_CDF_cons)
        fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f     '%m2logL + '\n')

print("CDF conservative done")


fresults.close()

print("***** scan finalized *****")

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
plt.xticks(np.linspace(-0.2, 0.4, 7))

data = np.genfromtxt(output)

## PDG
# Getting the data

x_PDG = data[0:9999,0]
y_PDG = data[0:9999,1]
z_PDG = data[0:9999,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2_PDG=[]
for z_el in z_PDG:
  z2_PDG.append(z_el-z_PDG.min())


# Interpolating the grid
xi_PDG = np.linspace(x_PDG.min(), x_PDG.max(), grid_subdivisions)
yi_PDG = np.linspace(y_PDG.min(), y_PDG.max(), grid_subdivisions)

X_PDG, Y_PDG = np.meshgrid(xi_PDG, yi_PDG)
Z_PDG = griddata((x_PDG, y_PDG), z2_PDG, (X_PDG, Y_PDG), method="linear")


# Plotting the 68% and 95% CL regions
cs1 = ax.contourf(xi_PDG,yi_PDG,Z_PDG,[10**(-10),2.3,5.99],colors=['#5d8aa8','#b0e0e6'], alpha=0.7)#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

## CDF
# Getting the data

x_CDF = data[10000:19999,0]
y_CDF = data[10000:19999,1]
z_CDF = data[10000:19999,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2_CDF=[]
for z_el in z_CDF:
  z2_CDF.append(z_el-z_CDF.min())


# Interpolating the grid
xi_CDF = np.linspace(x_CDF.min(), x_CDF.max(), grid_subdivisions)
yi_CDF = np.linspace(y_CDF.min(), y_CDF.max(), grid_subdivisions)

X_CDF, Y_CDF = np.meshgrid(xi_CDF, yi_CDF)
Z_CDF = griddata((x_CDF, y_CDF), z2_CDF, (X_CDF, Y_CDF), method="linear")


# Plotting the 68% and 95% CL regions
cs2 = ax.contourf(xi_CDF,yi_CDF,Z_CDF,[10**(-10),2.3,5.99],colors=['#e4717a','#ffb7c5'], alpha=0.7)#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

### CDF standard
## Getting the data

#x_CDF_std = data[20000:29999,0]
#y_CDF_std = data[20000:29999,1]
#z_CDF_std = data[20000:29999,2]


## Substracting the -2LogL minimum to form Delta(-2LogL)
#z2_CDF_std=[]
#for z_el in z_CDF_std:
#  z2_CDF_std.append(z_el-z_CDF_std.min())


## Interpolating the grid
#xi_CDF_std = np.linspace(x_CDF_std.min(), x_CDF_std.max(), grid_subdivisions)
#yi_CDF_std = np.linspace(y_CDF_std.min(), y_CDF_std.max(), grid_subdivisions)

#X_CDF_std, Y_CDF_std = np.meshgrid(xi_CDF_std, yi_CDF_std)
#Z_CDF_std = griddata((x_CDF_std, y_CDF_std), z2_CDF_std, (X_CDF_std, Y_CDF_std), method="linear")


## Plotting the 68% and 95% CL regions
#ax.contourf(xi_CDF_std,yi_CDF_std,Z_CDF_std,[10**(-10),2.3,5.99],colors=['#ff7369','#ff8c69'], alpha=0.5)#, \
#              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

### CDF conservative
## Getting the data

#x_CDF_cons = data[30000:39999,0]
#y_CDF_cons = data[30000:39999,1]
#z_CDF_cons = data[30000:39999,2]


## Substracting the -2LogL minimum to form Delta(-2LogL)
#z2_CDF_cons=[]
#for z_el in z_CDF_cons:
#  z2_CDF_cons.append(z_el-z_CDF_cons.min())


## Interpolating the grid
#xi_CDF_cons = np.linspace(x_CDF_cons.min(), x_CDF_cons.max(), grid_subdivisions)
#yi_CDF_cons = np.linspace(y_CDF_cons.min(), y_CDF_cons.max(), grid_subdivisions)

#X_CDF_cons, Y_CDF_cons = np.meshgrid(xi_CDF_cons, yi_CDF_cons)
#Z_CDF_cons = griddata((x_CDF_cons, y_CDF_cons), z2_CDF_cons, (X_CDF_cons, Y_CDF_cons), method="linear")


## Plotting the 68% and 95% CL regions
#ax.contourf(xi_CDF_cons,yi_CDF_cons,Z_CDF_cons,[10**(-10),2.3,5.99],colors=['#321fbc','#4c1fbc'], alpha=0.5)#, \
#              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])




#proxy1 = [plt.Rectangle((1, 1), 2, 2, fc=pc.get_facecolor()[0]) for pc in
#cs1.collections]

#proxy2 = [plt.Rectangle((1, 1), 2, 2, fc=pc.get_facecolor()[0]) for pc in
#cs2.collections]

#ax.legend([proxy1[0],proxy1[1], proxy2[0], proxy2[1]], ["1", "2", "3", "4"])



ax.set_aspect(1)

# Standard Model 
plt.plot([0], [0], '*',markersize=8, color = 'black', label="SM prediction")
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'S',fontsize=16)
plt.ylabel(r'T',fontsize=16)
plt.text(0.05, 0.270, r'CDF 2022', fontsize=11, rotation=40, color='#e4717a')
plt.text(-0.13, 0.015, r'PDG 2021', fontsize=11, rotation=47, color='#5d8aa8')
plt.text(-0.15, 0.45, r'fixed U = 0', fontsize=11)

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", validation_dir)
print("***** done *****")
