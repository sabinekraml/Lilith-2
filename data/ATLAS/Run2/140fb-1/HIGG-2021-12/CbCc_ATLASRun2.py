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

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "data/latestRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.38

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/cbcc-run2-HIGG-2021-12.out"
outputplot = "results/cbcc-run2-HIGG-2021-12.pdf"

# Scan ranges
CB_min = -3.75
CB_max = 3.75
CC_min = -25
CC_max = 20

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CB=1, CC=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CB, CC"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="tt">1</C>
  <C to="bb">%(CB)s</C>
  <C to="cc">%(CC)s</C>
  <C to="tautau">1</C>
  <C to="ZZ">1</C>
  <C to="WW">1</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CB':CB, 'CC':CC, 'precision':precision}
        
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


######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for CB in np.linspace(CB_min, CB_max, grid_subdivisions):
    fresults.write('\n')
    for CC in np.linspace(CC_min, CC_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, CB=CB, CC=CC, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CBmin = CB
            CCmin = CC
        fresults.write('%.5f    '%CB +'%.5f    '%CC + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CB, CC, -2logL_min = ", CBmin, CCmin, m2logLmin)

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=15, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)


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
dataload = open('data/ATLAS/Run2/140fb-1/HIGG-2021-12/officialData/HIGG-2021-12.csv','r')
dorix = []
doriy = []
for line in dataload:
  fdat = line.split(',')
  dorix.append(float(fdat[0]))
  doriy.append(float(fdat[1]))
 
# Plotting the 68%, 95% and 99.7% CL regions
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CB_max-CB_min)/(CC_max-CC_min))

# best fit point
plt.plot([CBmin],[CCmin], '*', c='w', ms=9)

# Standard Model 
plt.plot([1],[1], '+', c='k', ms=10, label = 'SM')

# best fit point - Official
plt.plot([-1.0180505415162449],[0.12244897959184442], 'P', c='k', ms=5, label = 'ATLAS official best fit')

# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=3,c='k',marker='o',label='ATLAS official')    
plt.legend(loc='lower right', scatterpoints = 3)

# Title, labels, color bar...
plt.title("  Lilith-2.1, ATLAS-HIGG-2021-12 validation" , fontsize=12, ha="center")
plt.xlabel(r'$C_B$',fontsize=20)
plt.ylabel(r'$C_C$',fontsize=20)
plt.text(-2.5, 16.5 , r'Exp. input type = vn', fontsize=12, ha = 'left')

fig.set_tight_layout(True)

#set aspect ratio to 0.618
ratio = 0.618
x_left, x_right = ax.get_xlim()
y_low, y_high = ax.get_ylim()
ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/results")
print("***** done *****")

