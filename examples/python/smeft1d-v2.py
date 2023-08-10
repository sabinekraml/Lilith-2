# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
import shutil
import re
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

print()    
print("".ljust(20,'*')+" STARTING ".ljust(35,"*"))
print()
lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith
lilithcalc = lilith.Lilith(verbose=False, timer=False)


print("     Import Lilith's core ".ljust(40,".")+" ok!".rjust(10,"."))

# ====================================================================== 
#			declaring input, scan region, output  
# ======================================================================

smeft_input = "data/smeft.list"  

#fit parameter 
fit_params = {}

fit_params['cHW'] = 'p1'

# Scan ranges
p1_min = -0.04
p1_max = 0.01

#number of grid steps in each of the two dimensions 
grid_subdivisions = 100


# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/SMEFT_1d.out"
outputplot = "results/SMEFT_1d.pdf"

officialplot = ''

# ======================================================================
# 			Loading and import input 
# ======================================================================

lilithcalc.readsmeftinput(smeft_input,fit_params)  
lilithcalc.readexpinput(lilithcalc.smeft_stxs)

# ======================================================================
#		    Scanning 
# ======================================================================

# Prepare output
fresults = open(output, 'w')

maxL = 10000
maxp1 = p1_min

print()
print("".ljust(15," ")+" Scanning ".rjust(15,"*")+"".rjust(5,"*"))
print()

for p1 in np.linspace(p1_min,p1_max,grid_subdivisions):
    fresults.write('\n')
    scan_value = {'p1':p1}
    lilithcalc.smeftmu_eval(scan_value)
    # ~ print(lilithcalc.smeft_user_mu)
    lilithcalc.computelikelihoodsmeft()
    result_l = lilithcalc.smeft_l
    if result_l < maxL:
        maxp1 = p1
        maxL = result_l 
    fresults.write('%.5f    '%p1 +'%.5f     '%result_l  +'\n')

fresults.close()
phrase = " at ("+list(fit_params.keys())[0]+") = "
print(" min(-2*logL) = ".rjust(35," "),"%10.5f"%maxL)
print(phrase.rjust(35," ")+" %10.5f"%maxp1)

lilithcalc.smeft_cleanlist()

#========================================================
#       plot
#========================================================

print()
print("".ljust(15," ")+" Plotting ".rjust(15,"*")+"".rjust(5,"*"))
print()

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
plt.tick_params(direction='in', labelsize=15, length=10, width=1.5)
plt.tick_params(which='minor', direction='in', length=6, width=1.0)


plt.ylim([-0.4,6])
plt.xlim([p1_min,p1_max])

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
# ~ y = data[:,1]
z = data[:,1]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())

plt.plot(x,z2,c="r",label='Lilith\'s fit ')

# Import Official data from file 
try:
    if officialplot=="":
        pass
    else:       
        dataload = open(officialplot,'r')
        dorix = []
        doriy = []
        for line in dataload:
            fdat = line.split(',')
            dorix.append(float(fdat[0]))
            doriy.append(float(fdat[1]))
        # Plotting the Offical contours 
        plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
        plt.legend(loc='lower right', scatterpoints = 3)  
except NameError: officialplot = None 

# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(p1_max+0.01*(p1_max - p1_min), 0.9, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(p1_max + 0.01*(p1_max - p1_min), 3.9, r'$2\sigma$', fontsize=10)
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

# Title, labels, color bar...
plt.title("  Lilith 2.2", fontsize=14, ha="left")
plt.xlabel(list(fit_params.keys())[0],fontsize=15)
plt.ylabel(r'$-2\ln\lambda$',fontsize=15)

# show plot or save plot
# ~ plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("     Results are stored in: <Lilith-directory>/results")


print()
print("".ljust(22,'-')+" DONE ".ljust(33,"-"))
print()
