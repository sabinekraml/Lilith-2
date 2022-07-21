import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CUCDCV/2d_simple/"

# Output files
output = validation_dir+"results/CVCD_140fb-1.out"
outputplot = validation_dir+"plots/CVCD_140fb-1.pdf"

######################################################################
# Plot routine
######################################################################
# Scan ranges
CD_min = 0.76
CD_max = 1.28
CV_min = 0.9
CV_max = 1.22

#Best Lilith fit point 
CVmin = 1.0422222222222222
CDmin = 1.0226262626262628
# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

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

ax.set_aspect((CV_max-CV_min)/(CD_max-CD_min))



# best Lilith fit point
plt.plot([CVmin],[CDmin], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
plt.title("CV-CD fit for 140 fb$^{-1}$ data", fontsize=12, ha="center")
plt.xlabel(r'$C_V$',fontsize=18)
plt.ylabel(r'$C_D$',fontsize=18)
#plt.text(0.82, 0.58, r'Data from CMS-HIG-19-015 Fig. 16 + Add. Fig. 3', fontsize=10)
#plt.text(0.84, 1.62, r'(Fig. 16 + Aux. Fig. 3)', fontsize=11)

plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", validation_dir)

