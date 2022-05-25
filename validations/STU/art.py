import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)


# Output files
#output = validation_dir+"mHmA_ST_mHpm.out"
#outputplot = validation_dir+"art/epeecosmique.pdf"
#output = validation_dir+"range/mHmA_ST_mHpm_60_200-700.out"
#outputplot = validation_dir+"art/araigneedemer.pdf"
output = validation_dir+"range/mHmA_ST_mHpm_30_200-700.out"
outputplot = validation_dir+"art/idkyet.pdf"


######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111, aspect='equal')


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
xi = np.linspace(x.min(), x.max(), np.sqrt(len(x)))
yi = np.linspace(y.min(), y.max(), np.sqrt(len(y)))

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")


# Plotting the 68% and 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),3.506,7.815],colors=['#ff3300','#ffa500',])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])


plt.axis('off')

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("results are stored in", lilith_dir + "validations/STU/art/")
print("***** done *****")
