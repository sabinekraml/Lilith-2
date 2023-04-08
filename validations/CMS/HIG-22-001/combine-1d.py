##################################################################
#
# Lilith routine example for (CV, CF)  plots
#
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Use the libraries matplotlib (plotting) and numpy (functions)
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "validations/CMS/HIG-22-001/LatestRun2-CMScombination-paperData.list" 
smpred_input = ""
smbin_corr_input = "" 

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "validations/CMS/HIG-22-001/output.out"

outputplot = "validations/CMS/HIG-22-001/CH-1d-combine.pdf"

# Scan ranges
CH_min = 0.92
CH_max = 1.05

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CV=1, CF=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CV, CF"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="tt">%(CF)s</C>
  <C to="bb">%(CF)s</C>
  <C to="cc">%(CF)s</C>
  <C to="tautau">%(CF)s</C>
  <C to="ZZ">%(CV)s</C>
  <C to="WW">%(CV)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CV':CV, 'CF':CF, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")


# Prepare output
fresults = open(output, 'w')

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False,timer=False)
# Read experimental data
lilithcalc.readexpinput(exp_input)
# Read SM prediction input and correlation 
lilithcalc.readsmpred(smpred_input)
lilithcalc.readsmcorr(smbin_corr_input)

######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for CH in np.linspace(CH_min, CH_max, grid_subdivisions):
    fresults.write('\n')
    myXML_user_input = usrXMLinput(hmass, CV=CH, CF=CH, precision=my_precision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    m2logL = lilithcalc.l
    if m2logL < m2logLmin:
        m2logLmin = m2logL
        CHmin = CH
    fresults.write('%.5f    '%CH + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CV, CF, -2logL_min = ", CHmin, m2logLmin)

#========================================================
#       plot
#========================================================
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
plt.tick_params(direction='out', labelsize=15, length=10, width=1.5)
plt.tick_params(which='minor', direction='out', length=6, width=1.0)

# Getting the data
data = np.genfromtxt(output)

x = data[:, 0]
y = data[:, 1]

# Substracting the -2LogL minimum to form Delta(-2LogL)
y2 = []
for y_el in y:
    y2.append(y_el - y.min())

# Plot the curve    
plt.plot(x,y2,c="r",label='Lilith\'s fit ')

# Import Official data from file 
#dataload = open('validations/ATLAS/HIGG-2021-23-Oct22/official-data/auxfig-1-total.csv','r')
#dorix = []
#doriy = []
#for line in dataload:
#    fdat = line.split(',')
#    dorix.append(float(fdat[0]))
#    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
#plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
#plt.legend(loc='lower right', scatterpoints = 3)  

# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(CH_max+0.1*(CH_max - CH_min), 0.8, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(CH_max + 0.1*(CH_max - CH_min), 3.8, r'$2\sigma$', fontsize=10)

# Title, labels, color bar...
plt.title("  Lilith-2.1, DB 22.x develop", fontsize=14, ha="left")
plt.xlabel(r'$cH\widetilde{W}B$',fontsize=15)
plt.ylabel(r'$-2\ln\lambda$',fontsize=15)
plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("***** done *****")
