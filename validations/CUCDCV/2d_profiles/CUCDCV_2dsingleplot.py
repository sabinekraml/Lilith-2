import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)

if (not os.path.exists("results")):
    os.mkdir("results")

if (not os.path.exists("plots")):
    os.mkdir("plots")


validation_dir = lilith_dir+"validations/CUCDCV/2d_profiles/"

# Output files
output = validation_dir+"results/CUCD-CV_2dprofiles_140fb-1_minuit.out"
outputplot = validation_dir+"plots/CUCD-CV_2dprofiles_140fb-1.pdf"





######################################################################
# Plot routine1
######################################################################
# Scan ranges
CU_min = 0.6
CU_max = 1.4
CD_min = 0.6
CD_max = 1.4

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

#Best Lilith fit point 
CUmin36 = 0.9636363636363636
CDmin36 = 0.9878787878787878

print("***** plotting data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(6,6) )
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

#ax.set_aspect((CU_max-CU_min)/(CD_max-CD_min))



# best Lilith fit point
plt.plot([CUmin36],[CDmin36], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
plt.title("CD-CU fit profiled over CV for 140 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_D$',fontsize=18)
plt.text(0.63, 0.63,"Lilith-2.1", fontsize=10)




plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", validation_dir)

