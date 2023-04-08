# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
#lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

#========================================================
#       plot 1 : aux-fig-12a cuH 
#========================================================

official_dat = "validations/ATLAS/HIGG-2018-28-SMEFT/official-data/auxfig-12a-cuH-even.csv"
output_ori = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12a-cuH-even.out"
output_no3 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12a-cuH-even-no3rd.out"
output_ln1 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12a-cuH-even-ln1.out"
output_ln2 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12a-cuH-even-ln2.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/1d-results/Lambda4-auxfig-12a-cuH.pdf"

print("***** Plotting 1/5 *****")

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

plt.ylim([-0.4,8])
plt.xlim([-50,35])

# Import Official data from file 
dataload = open(official_dat,'r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='m',marker='o',label='ATLAS official ')    
 
dataload.close()

# Full-terms data plot
data = np.genfromtxt(output_ori)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-',c='orange', linewidth=1, label='full terms')

# No 3rd-order-terms data plot
data = np.genfromtxt(output_no3)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='no 3rd-order-terms')

# Linear-only-data plot
data = np.genfromtxt(output_ln1)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-b', linewidth=1, label='only linear terms')

# non-linear-acceptance data plot
data = np.genfromtxt(output_ln2)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':g', linewidth=2, label='non-linear acceptance')



plt.legend(loc='upper right', scatterpoints = 3) 



# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(-45, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(-45, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith-2.2, data from ATLAS HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$c_{uH}$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("        .... done")

#========================================================
#       plot 2 : aux-fig-12a cuH 
#========================================================

official_dat = "validations/ATLAS/HIGG-2018-28-SMEFT/official-data/auxfig-12b-cHG-even.csv"
output_ori = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12b-cHG-even.out"
output_no3 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12b-cHG-even-no3rd.out"
output_ln1 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12b-cHG-even-ln1.out"
output_ln2 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12b-cHG-even-ln2.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/1d-results/Lambda4-auxfig-12b-cHG.pdf"

print("***** Plotting 2/5 *****")

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
plt.tick_params(direction='in', labelsize=10, length=7, width=1.0)
plt.tick_params(which='minor', direction='in', length=4, width=0.6)

plt.ylim([-0.4,8])
plt.xlim([-0.01,0.01])

# Import Official data from file 
dataload = open(official_dat,'r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
 
dataload.close()

# Full-terms data plot
data = np.genfromtxt(output_ori)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-',c='orange', linewidth=1, label='full terms')

# No 3rd-order-terms data plot
data = np.genfromtxt(output_no3)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='no 3rd-order-terms')

# Linear-only-data plot
data = np.genfromtxt(output_ln1)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-b', linewidth=1, label='only linear terms')

# non-linear-acceptance data plot
data = np.genfromtxt(output_ln2)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':g', linewidth=2, label='non-linear acceptance')



plt.legend(loc='upper right', scatterpoints = 3) 



# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(0.0085, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(0.0085, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith-2.2, data from ATLAS HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$c_{HG}$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("        .... done")

#========================================================
#       plot 3 : aux-fig-12c cHW 
#========================================================

official_dat = "validations/ATLAS/HIGG-2018-28-SMEFT/official-data/auxfig-12c-cHW-even.csv"
output_ori = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12c-cHW-even.out"
output_no3 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12c-cHW-even-no3rd.out"
output_ln1 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12c-cHW-even-ln1.out"
output_ln2 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12c-cHW-even-ln2.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/1d-results/Lambda4-auxfig-12c-cHW.pdf"

print("***** Plotting 3/5 *****")

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
plt.tick_params(direction='in', labelsize=10, length=7, width=1.0)
plt.tick_params(which='minor', direction='in', length=4, width=0.6)

plt.ylim([-0.4,8])
plt.xlim([-4.5,4.5])

# Import Official data from file 
dataload = open(official_dat,'r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
 
dataload.close()

# Full-terms data plot
data = np.genfromtxt(output_ori)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-',c='orange', linewidth=1, label='full terms')

# No 3rd-order-terms data plot
data = np.genfromtxt(output_no3)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='no 3rd-order-terms')

# Linear-only-data plot
data = np.genfromtxt(output_ln1)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-b', linewidth=1, label='only linear terms')

# non-linear-acceptance data plot
data = np.genfromtxt(output_ln2)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':g', linewidth=2, label='non-linear acceptance')



plt.legend(loc='upper right', scatterpoints = 3) 



# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(-4.3, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(-4.3, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith-2.2, data from ATLAS HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$c_{HW}$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("        .... done")


#========================================================
#       plot 4 : aux-fig-12d cHB 
#========================================================

official_dat = "validations/ATLAS/HIGG-2018-28-SMEFT/official-data/auxfig-12d-cHB-even.csv"
output_ori = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12d-cHB-even.out"
output_no3 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12d-cHB-even-no3rd.out"
output_ln1 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12d-cHB-even-ln1.out"
output_ln2 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12d-cHB-even-ln2.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/1d-results/Lambda4-auxfig-12d-cHB.pdf"

print("***** Plotting 4/5 *****")

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
plt.tick_params(direction='in', labelsize=10, length=7, width=1.0)
plt.tick_params(which='minor', direction='in', length=4, width=0.6)

plt.ylim([-0.4,8])
plt.xlim([-4,4])

# Import Official data from file 
dataload = open(official_dat,'r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
 
dataload.close()

# Full-terms data plot
data = np.genfromtxt(output_ori)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-',c='orange', linewidth=1, label='full terms')

# No 3rd-order-terms data plot
data = np.genfromtxt(output_no3)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='no 3rd-order-terms')

# Linear-only-data plot
data = np.genfromtxt(output_ln1)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-b', linewidth=1, label='only linear terms')

# non-linear-acceptance data plot
data = np.genfromtxt(output_ln2)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':g', linewidth=2, label='non-linear acceptance')



plt.legend(loc='upper right', scatterpoints = 3) 



# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(-3.8, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(-3.8, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith-2.2, data from ATLAS HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$c_{HB}$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("        .... done")


#========================================================
#       plot 5 : aux-fig-12c cHW 
#========================================================

official_dat = "validations/ATLAS/HIGG-2018-28-SMEFT/official-data/auxfig-12e-cHWB-even.csv"
output_ori = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12e-cHWB-even.out"
output_no3 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12e-cHWB-even-no3rd.out"
output_ln1 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12e-cHWB-even-ln1.out"
output_ln2 = "validations/ATLAS/HIGG-2018-28-SMEFT/output/auxfig-12e-cHWB-even-ln2.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/1d-results/Lambda4-auxfig-12e-cHWB.pdf"

print("***** Plotting 5/5 *****")

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
plt.tick_params(direction='in', labelsize=10, length=7, width=1.0)
plt.tick_params(which='minor', direction='in', length=4, width=0.6)

plt.ylim([-0.4,8])
plt.xlim([-4.5,6])

# Import Official data from file 
dataload = open(official_dat,'r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
 
dataload.close()

# Full-terms data plot
data = np.genfromtxt(output_ori)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-',c='orange', linewidth=1, label='full terms')

# No 3rd-order-terms data plot
data = np.genfromtxt(output_no3)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':r', linewidth=2, label='no 3rd-order-terms')

# Linear-only-data plot
data = np.genfromtxt(output_ln1)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,'-b', linewidth=1, label='only linear terms')

# non-linear-acceptance data plot
data = np.genfromtxt(output_ln2)
x = data[:, 0]
y = data[:, 1]	# no acceptance
z = data[:, 2]	# with acceptance
# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())
plt.plot(x,z2,':g', linewidth=2, label='non-linear acceptance')



plt.legend(loc='upper right', scatterpoints = 3) 



# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(-4.3, 1.1, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(-4.3, 4.1, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith-2.2, data from ATLAS HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$c_{HWB}$',fontsize=15)
plt.ylabel(r'$-2\ln\Lambda$',fontsize=15)
#plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
#plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("        .... done")
