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
fit_params['cHB'] = 'p2'

# Scan ranges
fit_range = {}

p1_min = -0.2
p1_max = 0.3
p2_min = -0.10
p2_max = 0.07

#number of grid steps in each of the two dimensions 
grid_subdivisions = 100

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/SMEFT_2d.out"
outputplot = "results/SMEFT_2d.pdf"

officialplot = ''

# for plotting axes label 
fit_params_list = list(fit_params.keys()) 
# ~ fit_params_list = list(fit_params.values()) 

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
maxp2 = p2_min

print()
print("".ljust(15," ")+" Scanning ".rjust(15,"*")+"".rjust(5,"*"))
print()

# Scan

for y in np.linspace(p2_min,p2_max,grid_subdivisions):
    fresults.write('\n')
    for x in np.linspace(p1_min,p1_max,grid_subdivisions):
        scan_value = {"p1":x, "p2":y}
        lilithcalc.smeftmu_eval(scan_value)
        lilithcalc.computelikelihoodsmeft()
        result_l = lilithcalc.smeft_l
        if result_l < maxL:
            maxp1 = x
            maxp2 = y
            maxL = result_l 
        fresults.write('%.5f    '%x +'%.5f    '%y +'%.5f     '%result_l  +'\n')

fresults.close()

phrase = " at ("+fit_params_list[0]+","+fit_params_list[1]+") = "
print(" min(-2*logL) = ".rjust(35," "),maxL)
print(phrase.rjust(35," ")+"( %.4f"%maxp1+" , %.4f"%maxp2+" )")

lilithcalc.smeft_cleanlist()

#========================================================
#       plot
#========================================================

print()
print("".ljust(15," ")+" Plotting ".rjust(15,"*")+"".rjust(5,"*"))
print()

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
plt.tick_params(direction='in', labelsize=12, length=10, width=1.5)
plt.tick_params(which='minor', direction='in', length=7, width=1.0)

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())
t2=[]  

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")

# Plotting the 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((p1_max-p1_min)/(p2_max-p2_min))

# best fit point
plt.plot([maxp1],[maxp2], '*', c='k', ms=5, zorder=1, label='Lilith best-fit')

# Standard Model 
plt.plot([0],[0], '+', c='k', ms=10, label='SM prediction')

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
        plt.scatter(dorix,doriy,s=3,c='k',marker='o',label='ATLAS official ',zorder=1)    
        plt.legend(loc='lower left', scatterpoints = 3) 
except NameError: officialplot = None 



# Title, labels, color bar...
plt.title("  Lilith-2.2", fontsize=8, ha="left")
plt.xlabel(fit_params_list[0],fontsize=15)
plt.ylabel(fit_params_list[1],fontsize=15)
plt.legend(loc='upper right')

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("     Results are stored in: <Lilith-directory>/results")


print()
print("".ljust(22,'-')+" DONE ".ljust(33,"-"))
print()
