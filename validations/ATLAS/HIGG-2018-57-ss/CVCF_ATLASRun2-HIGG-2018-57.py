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
exp_input = "validations/ATLAS/HIGG-2018-57-ss/LatestRun2.list"
smpred_input = "validations/ATLAS/HIGG-2018-57-ss/SM-prediction-ss.txt"
smbin_corr_input = "" 

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.38

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "validations/ATLAS/HIGG-2018-57-ss/CVCF-HIGG-2018-57.out"
#outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-noSMerr.pdf"
#outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-SMerr.pdf"
#outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-cov-mod.pdf"
#outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-approx-g1.pdf"
#outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-approx-g2.pdf"
outputplot = "validations/ATLAS/HIGG-2018-57-ss/figs/CVCF-HIGG-2018-57-ss-wSMprediction.pdf"

# Scan ranges
CV_min = 0.9
CV_max = 1.17
CF_min = 0.75
CF_max = 1.45
# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 50

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
lilithcalc.readsmpred(smpred_input)
lilithcalc.readsmcorr(smbin_corr_input)
######################################################################
# Scan routine
######################################################################

m2logLmin=10000000000
max=-1

print("***** running scan *****")

for CV in np.linspace(CV_min, CV_max, grid_subdivisions):
    fresults.write('\n')
    for CF in np.linspace(CF_min, CF_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, CV=CV, CF=CF, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CVmin = CV
            CFmin = CF
        fresults.write('%.5f    '%CV +'%.5f    '%CF + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CV, CF, -2logL_min = ", CVmin, CFmin, m2logLmin)

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 5
matplotlib.rcParams['ytick.major.pad'] = 5

fig = plt.figure()
ax = fig.add_subplot(111)

ax.locator_params(axis='y', nbins=10)
ax.locator_params(axis='x', nbins=10)

#ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

#ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=1.2)
plt.tick_params(which='minor', direction='in', length=5, width=0.7)



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

# Import Official data from file 
dataload = open('validations/ATLAS/HIGG-2018-57-ss/HIGG-2018-57-official.csv','r')
dorix = []
doriy = []
for line in dataload:
  fdat = line.split(',')
  dorix.append(float(fdat[0]))
  doriy.append(float(fdat[1]))
   
# Plotting the 68%, 95% and 99.7% CL regions
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CV_max-CV_min)/(CF_max-CF_min))

# best fit point
plt.plot([CVmin],[CFmin], '*', c='w', ms=10)

# Standard Model 
plt.plot([1],[1], '+', c='k', ms=10)

# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=3,c='b',marker='o',label='ATLAS official ')    
plt.legend(loc='lower right', scatterpoints = 3) 

# best fit point
#plt.plot([1.053485254691689],[1.0492700729927007], 'o', c='b', ms=3)

# Title, labels, color bar...
plt.title("  Lilith 2.2 - ATLAS HIGG-2018-57" , fontsize=18, ha="center")
plt.xlabel(r'$C_V$',fontsize=18)
plt.ylabel(r'$C_F$',fontsize=18)

plt.text(0.915, 1.41, r'Exp. input type = SS', fontsize=12, ha = 'left')
#plt.text(0.915, 1.37, r'include theoretical error', fontsize=12, ha = 'left')
#plt.text(0.915, 1.33, r'correlation approx, gamma = 2.0', fontsize=12, ha = 'left')
plt.text(0.915, 1.37, r'include theoretical error', fontsize=12, ha = 'left')
plt.text(0.915, 1.33, r'add Theo. err. from twiki.cern.ch', fontsize=12, ha = 'left')
plt.text(0.915, 1.29, r'distribution = vn', fontsize=12, ha = 'left')

fig.set_tight_layout(True)

#set aspect ratio to 0.618
#ratio = 0.85
#x_left, x_right = ax.get_xlim()
#y_low, y_high = ax.get_ylim()
#ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot,bbox_inches='tight')

print("results are stored in", lilith_dir + "/results")
print("***** done *****")

