import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 5
matplotlib.rcParams['ytick.major.pad'] = 5

fig = plt.figure()
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=12, length=10, width=1.5)
plt.tick_params(which='minor', direction='in', length=7, width=1.0)



# Getting the data
data = np.genfromtxt("Mathpoints.txt")

x = data[:,0]
y = data[:,1]
z = data[:,2]
#t = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())
t2=[]  
#for t_el in t:
#  t2.append(t_el-t.min())

# Interpolating the grid
grid_subdivisions = 100
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")
#T = griddata((x, y), t2, (X, Y), method="linear")

# Plotting the 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),5.99],colors=['#ffa500'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],zorder=-1, alpha=1)
#ax.contourf(xi,yi,T,[10**(-10),5.99],colors=['#bd8fbe'], \
#              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],zorder=-0.5, alpha=0.7)	      

#ax.set_aspect((cx_max-cx_min)/(cy_max-cy_min))

# best fit point
#plt.plot([maxpx_ac],[maxpy_ac], '*', c='r', ms=10,zorder=0)
#plt.plot([maxpx],[maxpy], '*', c='b', ms=10,zorder=0)

# Standard Model 
plt.plot([0],[0], '+', c='k', ms=10)

# Import Official data from file 
dataload = open('official-data/fig18a-cHW-cHB-odd.csv','r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=3,c='k',marker='o',label='ATLAS official ',zorder=1)    
plt.legend(loc='lower left', scatterpoints = 3)  



# Title, labels, color bar...
plt.title("  Lilith-2.1, DB 22.x develop", fontsize=8, ha="left")
plt.xlabel(r'$cH\widetilde{W}$',fontsize=15)
plt.ylabel(r'$cH\widetilde{B}$',fontsize=15)
plt.text(-0, 3.5, r'Data from ATLAS HIGG-2018-28', fontsize=8)
plt.text(-0, 3.2, r'$H\to ZZ^*\to 4\ell$', fontsize=8)

fig.set_tight_layout(True)

plt.show()

# Saving figure (.pdf)
#fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/validations/ATLAS/HIGG-2018-28-SMEFT/2d-results")
print("***** done *****")
