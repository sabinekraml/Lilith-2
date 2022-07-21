import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))+"/"
sys.path.append(lilith_dir)

validation_dir = lilith_dir+"validations/CUCDCV/2d_profiles/"


# Output files
output36_1 = validation_dir+"results/CDCV-CU_2dprofiles_36fb-1_minuit.out"
output140_1 = validation_dir+"results/CDCV-CU_2dprofiles_140fb-1_minuit.out"
output36_2 = validation_dir+"results/CUCD-CV_2dprofiles_36fb-1_minuit.out"
output140_2 = validation_dir+"results/CUCD-CV_2dprofiles_140fb-1_minuit.out"
output36_3 = validation_dir+"results/CUCV-CD_2dprofiles_36fb-1_minuit.out"
output140_3 = validation_dir+"results/CUCV-CD_2dprofiles_140fb-1_minuit.out"

outputplot = "profilingfits_colormap.pdf"





######################################################################
# Plot routine1
######################################################################
# Scan ranges
CD_min1 = 0.7
CD_max1 = 1.4
CV_min1 = 0.8
CV_max1 = 1.3

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

#Best Lilith fit point 
CDmin36_1 = 1.0606060606060606
CVmin36_1 = 1.077777777777778

print("***** plotting 36fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(48,12) )
ax = fig.add_subplot(231)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_36_1 = np.genfromtxt(output36_1)

x36_1 = data_36_1[:,0]
y36_1 = data_36_1[:,1]
z36_1 = data_36_1[:,2]
t36_1 = data_36_1[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z236_1=[]
for z_el in z36_1:
  z236_1.append(z_el-z36_1.min())


# Interpolating the grid
xi36_1 = np.linspace(x36_1.min(), x36_1.max(), grid_subdivisions)
yi36_1 = np.linspace(y36_1.min(), y36_1.max(), grid_subdivisions)

X36_1, Y36_1 = np.meshgrid(xi36_1, yi36_1)
#print("X36 après meshgrid:", X36_1)
Z36_1 = griddata((x36_1, y36_1), z236_1, (X36_1, Y36_1), method="linear")
#print("griddata:", Z36_1)

# Plotting the 68% and 95% CL regions
ax.contour(xi36_1,yi36_1,Z36_1,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CD_max1-CD_min1)/(CV_max1-CV_min1))

sc = ax.scatter(x36_1, y36_1, c=t36_1, cmap="jet_r")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$C_U$", fontsize=18)

# best Lilith fit point
plt.plot([CDmin36_1],[CVmin36_1], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CD-CV fit profiled over CU for 36 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_D$',fontsize=18)
plt.ylabel(r'$C_V$',fontsize=18)
plt.text(0.64, 0.83,"Lilith-2.1", fontsize=10)



######################################################################
# Plot routine2
######################################################################

#Best Lilith fit point 
CDmin140_1 = 0.9878787878787878
CVmin140_1 = 1.0373737373737375

print("***** plotting 140_1fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

ax = fig.add_subplot(234)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_140_1 = np.genfromtxt(output140_1)

x140_1 = data_140_1[:,0]
y140_1 = data_140_1[:,1]
z140_1 = data_140_1[:,2]
t140_1 = data_140_1[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2140_1=[]
for z_el in z140_1:
  z2140_1.append(z_el-z140_1.min())


# Interpolating the grid
xi140_1 = np.linspace(x140_1.min(), x140_1.max(), grid_subdivisions)
yi140_1 = np.linspace(y140_1.min(), y140_1.max(), grid_subdivisions)

X140_1, Y140_1 = np.meshgrid(xi140_1, yi140_1)
Z140_1 = griddata((x140_1, y140_1), z2140_1, (X140_1, Y140_1), method="linear")


# Plotting the 68% and 95% CL regions
ax.contour(xi140_1,yi140_1,Z140_1,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

sc2 = ax.scatter(x140_1, y140_1, c=t140_1, cmap="jet_r")
cbar = fig.colorbar(sc2,fraction=0.046, pad=0.04)
cbar.set_label("$C_U$", fontsize=15)

ax.set_aspect((CD_max1-CD_min1)/(CV_max1-CV_min1))



# best Lilith fit point
plt.plot([CDmin140_1],[CVmin140_1], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CD-CV fit profiled over CU for 140 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_D$',fontsize=18)
plt.ylabel(r'$C_V$',fontsize=18)
plt.text(0.64, 0.83,"Lilith-2.1", fontsize=10)



######################################################################
# Plot routine3
######################################################################
# Scan ranges
CU_min2 = 0.6
CU_max2 = 1.4
CD_min2 = 0.6
CD_max2 = 1.5


#Best Lilith fit point 
CUmin36_2 = 1.0444444444444443
CDmin36_2 = 1.0606060606060606

print("***** plotting 36fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8


ax = fig.add_subplot(232)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_36_2 = np.genfromtxt(output36_2)

x36_2 = data_36_2[:,0]
y36_2 = data_36_2[:,1]
z36_2 = data_36_2[:,2]
t36_2 = data_36_2[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z236_2=[]
for z_el in z36_2:
  z236_2.append(z_el-z36_2.min())


# Interpolating the grid
xi36_2 = np.linspace(x36_2.min(), x36_2.max(), grid_subdivisions)
yi36_2 = np.linspace(y36_2.min(), y36_2.max(), grid_subdivisions)

X36_2, Y36_2 = np.meshgrid(xi36_2, yi36_2)
Z36_2 = griddata((x36_2, y36_2), z236_2, (X36_2, Y36_2), method="linear")


# Plotting the 68% and 95% CL regions
ax.contour(xi36_2,yi36_2,Z36_2,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max2-CU_min2)/(CD_max2-CD_min2))


sc3 = ax.scatter(x36_2, y36_2, c=t36_2, cmap="jet_r")
cbar = fig.colorbar(sc3,fraction=0.046, pad=0.04)
cbar.set_label("$C_V$", fontsize=15)


# best Lilith fit point
plt.plot([CUmin36_2],[CDmin36_2], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CD-CU fit profiled over CV for 36 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_D$',fontsize=18)
plt.text(0.63, 0.63,"Lilith-2.1", fontsize=10)



######################################################################
# Plot routine4
######################################################################

#Best Lilith fit point 
CUmin140_2 = 0.9636363636363636
CDmin140_2 = 0.9878787878787878

print("***** plotting 140_2fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

ax = fig.add_subplot(235)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_140_2 = np.genfromtxt(output140_2)

x140_2 = data_140_2[:,0]
y140_2 = data_140_2[:,1]
z140_2 = data_140_2[:,2]
t140_2 = data_140_2[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2140_2=[]
for z_el in z140_2:
  z2140_2.append(z_el-z140_2.min())


# Interpolating the grid
xi140_2 = np.linspace(x140_2.min(), x140_2.max(), grid_subdivisions)
yi140_2 = np.linspace(y140_2.min(), y140_2.max(), grid_subdivisions)

X140_2, Y140_2 = np.meshgrid(xi140_2, yi140_2)
Z140_2 = griddata((x140_2, y140_2), z2140_2, (X140_2, Y140_2), method="linear")


# Plotting the 68% and 95% CL regions
ax.contour(xi140_2,yi140_2,Z140_2,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max2-CU_min2)/(CD_max2-CD_min2))

sc4 = ax.scatter(x140_2, y140_2, c=t140_2, cmap="jet_r")
cbar = fig.colorbar(sc4,fraction=0.046, pad=0.04)
cbar.set_label("$C_V$", fontsize=15)

# best Lilith fit point
plt.plot([CUmin140_2],[CDmin140_2], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CD-CU fit profiled over CV for 140 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_D$',fontsize=18)
plt.text(0.63, 0.63,"Lilith-2.1", fontsize=10)




######################################################################
# Plot routine5
######################################################################

# Scan ranges
CU_min3 = 0.75
CU_max3 = 1.25
CV_min3 = 0.8
CV_max3 = 1.3

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

#Best Lilith fit point 
CUmin36_3 = 1.0429292929292928
CVmin36_3 = 1.077777777777778


print("***** plotting 36fb-1_3 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

ax = fig.add_subplot(233)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_36_3 = np.genfromtxt(output36_3)

x36_3 = data_36_3[:,0]
y36_3 = data_36_3[:,1]
z36_3 = data_36_3[:,2]
t36_3 = data_36_3[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z236_3=[]
for z_el in z36_3:
  z236_3.append(z_el-z36_3.min())


# Interpolating the grid
xi36_3 = np.linspace(x36_3.min(), x36_3.max(), grid_subdivisions)
yi36_3 = np.linspace(y36_3.min(), y36_3.max(), grid_subdivisions)

X36_3, Y36_3 = np.meshgrid(xi36_3, yi36_3)
Z36_3 = griddata((x36_3, y36_3), z236_3, (X36_3, Y36_3), method="linear")


# Plotting the 68% and 95% CL regions
ax.contour(xi36_3,yi36_3,Z36_3,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max3-CU_min3)/(CV_max3-CV_min3))

sc5 = ax.scatter(x36_3, y36_3, c=t36_3, cmap="jet_r")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$C_D$", fontsize=15)


# best Lilith fit point
plt.plot([CUmin36_3],[CVmin36_3], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CU-CV fit profiled over CD for 36 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_V$',fontsize=18)
plt.text(0.78, 0.825,"Lilith-2.1", fontsize=10)


######################################################################
# Plot routine6
######################################################################

#Best Lilith fit point 
CUmin140_3 = 0.9621212121212122
CVmin140_3 = 1.0373737373737375

print("***** plotting 140_3fb-1 data fit *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

ax = fig.add_subplot(236)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data_140_3 = np.genfromtxt(output140_3)

x140_3 = data_140_3[:,0]
y140_3 = data_140_3[:,1]
z140_3 = data_140_3[:,2]
t140_3 = data_140_3[:,3]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2140_3=[]
for z_el in z140_3:
  z2140_3.append(z_el-z140_3.min())


# Interpolating the grid
xi140_3 = np.linspace(x140_3.min(), x140_3.max(), grid_subdivisions)
yi140_3 = np.linspace(y140_3.min(), y140_3.max(), grid_subdivisions)

X140_3, Y140_3 = np.meshgrid(xi140_3, yi140_3)
Z140_3 = griddata((x140_3, y140_3), z2140_3, (X140_3, Y140_3), method="linear")


# Plotting the 68% and 95% CL regions
ax.contour(xi140_3,yi140_3,Z140_3,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CU_max3-CU_min3)/(CV_max3-CV_min3))

sc6 = ax.scatter(x140_3, y140_3, c=t140_3, cmap="jet_r")
cbar = fig.colorbar(sc6,fraction=0.046, pad=0.04)
cbar.set_label("$C_D$", fontsize=15)

# best Lilith fit point
plt.plot([CUmin140_3],[CVmin140_3], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
#plt.title("CU-CV fit profiled over CD for 140 fb$^{-1}$ data", fontsize=13, ha="center")
plt.xlabel(r'$C_U$',fontsize=18)
plt.ylabel(r'$C_V$',fontsize=18)
plt.text(0.78, 0.825,"Lilith-2.1", fontsize=10)




#fig.tight_layout()    
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.65, hspace=0.2)          
plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", outputplot)

