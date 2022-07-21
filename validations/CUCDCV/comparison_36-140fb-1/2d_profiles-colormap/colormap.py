import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))+"/"
sys.path.append(lilith_dir)


validation_dir = lilith_dir+"validations/CUCDCV/"

# Output files

output140 = validation_dir+"CUCV-CD_2dprofiles_140fb-1zoom_minuit.out"
outputplot = validation_dir+"CUCV-CD_colormap_zoom.pdf"





######################################################################
# Plot routine
######################################################################
# Scan ranges
CU_min = 0.75
CU_max = 1.15
CV_min = 0.9
CV_max = 1.2

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

#Best Lilith fit point 
CUmin36 = 0.9641414141414141 
CDmin36 = 1.0388888888888888


print("***** plotting data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(8,8) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(output140)

x36 = data[:,0]
y36 = data[:,1]
z36 = data[:,2]
t36 = data[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z236=[]
for z_el in z36:
  z236.append(z_el-z36.min())


# Interpolating the grid
xi36 = np.linspace(x36.min(), x36.max(), grid_subdivisions)
yi36 = np.linspace(y36.min(), y36.max(), grid_subdivisions)

X36, Y36 = np.meshgrid(xi36, yi36)
Z36 = griddata((x36, y36), z236, (X36, Y36), method="linear")


# Plotting the 68% and 95% CL regions
cs=ax.contour(xi36,yi36,Z36,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'], linestyles=["dashed","dashed"])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max-CU_min)/(CV_max-CV_min))
#ax.set_xlim([0.75,1.15])
#ax.set_ylim([0.9,1.15])


# best Lilith fit point
plt.plot([CUmin36],[CDmin36], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1],'+',markersize=9, color = '#252850', label="SM prediction")



# Title, labels, color bar...
plt.title("CU-CV fit profiled over CD for 140 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_V$',fontsize=18)
plt.text(1.16, 0.82,"Lilith-2.1", fontsize=10)

sc = ax.scatter(x36, y36, c=t36, cmap="Spectral")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$C_D$", fontsize=15)
#labels = ["68%", "95%"] 
#for i in range(len(labels)): 
#    cs.collections[i].set_label(labels[i])
    
plt.legend(loc='upper left')   
plt.savefig(outputplot)
plt.show()
