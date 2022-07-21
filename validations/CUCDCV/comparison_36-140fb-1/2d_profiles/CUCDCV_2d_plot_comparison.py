import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))+"/"
sys.path.append(lilith_dir)


validation_dir = lilith_dir+"validations/CUCDCV/2d_profiles/"

if (not os.path.exists("plots")):
    os.mkdir("plots")
    
# Output files
output36 = validation_dir+"results/CUCD-CV_2dprofiles_36fb-1_minuit.out"
output140 = validation_dir+"results/CUCD-CV_2dprofiles_140fb-1_minuit.out"
outputplot = "plots/CUCD-CV_2dprofiles.pdf"





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
CUmin36 = 1.0444444444444443
CDmin36 = 1.0606060606060606

print("***** plotting 36fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(12,8) )
ax = fig.add_subplot(121)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(output36)

x36 = data[:,0]
y36 = data[:,1]
z36 = data[:,2]


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
ax.contourf(xi36,yi36,Z36,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max-CU_min)/(CD_max-CD_min))



# best Lilith fit point
plt.plot([CUmin36],[CDmin36], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
plt.title("CD-CU fit profiled over CV for 36 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_D$',fontsize=18)
plt.text(0.63, 0.63,"Lilith-2.1", fontsize=10)



######################################################################
# Plot routine2
######################################################################

#Best Lilith fit point 
CUmin140 = 0.9636363636363636
CDmin140 = 0.9878787878787878

print("***** plotting 140fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

ax = fig.add_subplot(122)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(output140)

x140 = data[:,0]
y140 = data[:,1]
z140 = data[:,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2140=[]
for z_el in z140:
  z2140.append(z_el-z140.min())


# Interpolating the grid
xi140 = np.linspace(x140.min(), x140.max(), grid_subdivisions)
yi140 = np.linspace(y140.min(), y140.max(), grid_subdivisions)

X140, Y140 = np.meshgrid(xi140, yi140)
Z140 = griddata((x140, y140), z2140, (X140, Y140), method="linear")


# Plotting the 68% and 95% CL regions
ax.contourf(xi140,yi140,Z140,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max-CU_min)/(CD_max-CD_min))



# best Lilith fit point
plt.plot([CUmin140],[CDmin140], '*', markersize=8, color = 'black', label = 'Lilith best fit')


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
print("results are stored in", outputplot)

