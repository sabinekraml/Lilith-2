##################################################################
#
# Plot routine for 2d likelihood on the $m_H$, $m_A$ with $m_Hpm$ minimized at each point
# Each point is in the 2-sigma interval, color code indicating the $mH_pm$ value 
# 
##################################################################

import sys, os
from scipy.interpolate import griddata
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/rangeminimize/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Scan ranges
mA_min = 200
mA_max = 2000
mH_min = 200
mH_max = 2000
mHpm_min = 200
mHpm_max = 2000

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

# Output files
output = validation_dir+"mHmA_ST_mHpm_100_200-2000.out"
outputplot = validation_dir+"mHmA_ST_mHpm_100_200-2000_colormHpm.pdf"


######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 4
matplotlib.rcParams['ytick.major.pad'] = 4

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=10, length=8, width=1.5)
#plt.minorticks_off()
plt.tick_params(which='minor', direction='in', length=5, width=0.5)


# Getting the data
data = np.genfromtxt(output)

x = data[0:-1,0]
y = data[0:-1,1]
z = data[0:-1,2]
zz = data[0:-1,3]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")

# Keeping only the 2-sigma points
j=-1
toremove = []
for i in z2:
    j+=1
    if i>7.815:
        toremove.append(j)

x = np.delete(x, toremove)
y = np.delete(y, toremove)
z2 = np.delete(z2, toremove)
zz = np.delete(zz, toremove)

# Plotting
sc = ax.scatter(x, y, c=zz, cmap="jet_r")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$m_{{H^{{\pm}}}}$[GeV]", fontsize=10)

ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

# best fit point
plt.plot([data[-1,0]],[data[-1,1]], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$m_H$[GeV]',fontsize=12)
plt.ylabel(r'$m_A$[GeV]',fontsize=12)
plt.text(mH_min + 100, mA_min + 350, r'Scatter plot in the $m_H$, $m_A$ plane with $m_{H^{\pm}}$ minimized at each point', fontsize=6)
plt.text(mH_min + 100, mA_min + 250, fr"Range of $m_{{H^{{\pm}}}}$ = ({mHpm_min},{mHpm_max}), with $\cos(\beta - \alpha) = 0$", fontsize=6)
plt.text(mH_min + 100, mA_min + 150, fr"Best point ($m_H, m_A$) = ({data[-1,0]:.0f}, {data[-1,1]:.0f}) with $\chi^{2}$ = {data[-1,2]:.3f} and $m_{{H^{{\pm}}}}$ = {data[-1,3]:.0f}", fontsize=6)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "validations/STU/")
print("***** done *****")
