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
from mpl_toolkits.mplot3d import Axes3D

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
STcorrelation = 0.9
SUcorrelation = -0.59
TUcorrelation = -0.85

# Output files
#output = validation_dir+"STU0PDG.out"
output = validation_dir+"STU.out"
outputplot = validation_dir+"STU.pdf"

# Scan ranges
S_min = -0.5
S_max = 0.5
T_min = -0.5
T_max = 0.5
U_min = -0.5
U_max = 0.5

CEN = np.array([Scen, Tcen, Ucen])
SIG = np.diag([Ssigma, Tsigma, Usigma])
COR = np.array(([1, STcorrelation, SUcorrelation],
                [STcorrelation, 1 , TUcorrelation],
                [SUcorrelation, TUcorrelation, 1]))
C = SIG.dot(COR).dot(SIG)
C_inv = np.linalg.inv(C)

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 5

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

def funcmatrix(X, CEN, SIG, COR):
    return C_inv.dot(X-CEN).dot((X-CEN).T)

######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for S in np.linspace(S_min, S_max, grid_subdivisions):
    fresults.write('\n')
    for T in np.linspace(T_min, T_max, grid_subdivisions):
        m2logL = func([S,T], Scen, Tcen, Ssigma, Ssigma, Tsigma, Tsigma, STcorrelation)
        for U in np.linspace(U_min, U_max, grid_subdivisions):
            m2logL = funcmatrix([S,T,U], CEN, SIG, COR)
            if m2logL < m2logLmin:
                m2logLmin = m2logL
                Smin = S
                Tmin = T
                Umin = U
            fresults.write('%.5f    '%S +'%.5f    '%T + '%.5f    '%U + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at S, T, U, -2logL_min = ", Smin, Tmin, Umin, m2logLmin)


#######################################################################
## Plot routine
#######################################################################

print("***** plotting *****")

# Preparing plot
#matplotlib.rcParams['xtick.major.pad'] = 8
#matplotlib.rcParams['ytick.major.pad'] = 8
#matplotlib.rcParams['ztick.major.pad'] = 8

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.yaxis.set_ticks_position('both')
#ax.yaxis.set_label_position('left')

#ax.xaxis.set_ticks_position('both')
#ax.xaxis.set_label_position('bottom')

#plt.minorticks_on()
#plt.tick_params(direction='in', labelsize=14, length=10, width=2)
#plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]
l = data[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
l2=[]
for l_el in l:
  l2.append(l_el-l.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)
zi = np.linspace(z.min(), z.max(), grid_subdivisions)

X, Y, Z = np.meshgrid(xi, yi, zi)
L = griddata((x, y, z), l2, (X, Y, Z), method="linear")

# Plotting the 68% and 95% CL regions
#ax.contourf(xi,yi,zi,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
#              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
ax.contourf(X,Y,Z,c='L2')
#ax.scatter(x,y,z,c=l2)

print("test")

# best fit point
plt.plot([Smin],[Tmin],[Umin], '+', markersize=8, color = 'white', label = 'CDF best fit')

# Standard Model 
plt.plot([0], [0], [0], '*',markersize=8, color = 'black', label="SM prediction")
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'S',fontsize=16)
plt.ylabel(r'T',fontsize=16)
ax.set_zlabel(r'U',fontsize=16)
#plt.text(-0.15, 0.45, r'fixed U = 0', fontsize=11)

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "validations/STU/")
print("***** done *****")
