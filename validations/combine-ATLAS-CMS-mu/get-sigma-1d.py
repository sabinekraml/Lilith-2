# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

official_dat = ''
output1 = "validations/combine-ATLAS-CMS-mu/output/combine-mu.out"
output2 = "validations/combine-ATLAS-CMS-mu/output/combine-CH.out"
outputplot = "validations/combine-ATLAS-CMS-mu/output/combine-graph.pdf"
print("***** Plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=10, length=7, width=1.0)
plt.tick_params(which='minor', direction='in', length=4, width=0.6)

plt.ylim([-0.4,6])
#plt.xlim([0.85,1.25])   #ATLAS
plt.xlim([0.9,1.15])

# Import Official data from file 
#dataload = open(official_dat,'r')
#dorix = []
#doriy = []
#for line in dataload:
#    fdat = line.split(',')
#    dorix.append(float(fdat[0]))
#    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
#plt.scatter(dorix,doriy,s=4,c='m',marker='o',label='ATLAS official ')    
#dataload.close()

# Plot output 1
print(" data 1: from Cmu ")
data = np.genfromtxt(output1)
x = data[:, 0]
z = data[:, 1]	
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
#take best-fit, left/right uncertainties    
min_value = min(z2)
min_index = z2.index(min_value)
z3 = []
for z_index in range(len(z2)):
    if z2[z_index] <= 1.0:
        z3.append(z_index)
neg_index = z3[0]
pos_index = z3[len(z3)-1] 
print("           best fit at          ",x[min_index])
print("           left uncertainty  = -",x[min_index] - x[neg_index])
print("           right uncertainty = +",x[pos_index] - x[min_index])  
plt.plot(x,z2,'-',c='orange', linewidth=1, label='input = Cmu')



# Plot output 2
print(" data 2: from C**2 ")
data = np.genfromtxt(output2)
x = data[:, 0]
z = data[:, 1]	
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='input = CH^2')
plt.legend(loc='upper center', scatterpoints = 3) 
#take best-fit, left/right uncertainties    
min_value = min(z2)
min_index = z2.index(min_value)
z3 = []
for z_index in range(len(z2)):
    if z2[z_index] <= 1.0:
        z3.append(z_index)
neg_index = z3[0]
pos_index = z3[len(z3)-1] 
print("           best fit at          ",x[min_index])
print("           left uncertainty  = -",x[min_index] - x[neg_index])
print("           right uncertainty = +",x[pos_index] - x[min_index])  


# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(1.07, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(1.07, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
#plt.title("  Lilith-2.2, data from ATLAS HIGG-2021-23", fontsize=14, ha="center")
plt.title("  Lilith-2.2, data from CMS HIG-22-001", fontsize=14, ha="center")
plt.xlabel(r'$mu$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(, , r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 1.22, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print(" Plotting ... done!")
