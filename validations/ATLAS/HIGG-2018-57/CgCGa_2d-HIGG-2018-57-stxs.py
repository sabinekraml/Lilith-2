###############################################################
#
# Lilith routine example for (C_gluon, C_gamma) plots
#
# To put in Lilith-2.X/examples/python/ folder 
# To run from /Lilith-2.X root folder
#
# Use the libraries matplotlib (plotting) and numpy (functions)
#
###############################################################

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
exp_input = "validations/ATLAS/HIGG-2018-57/LatestRun2.list"
# Sm predictions     
smpred_input = "validations/ATLAS/HIGG-2018-57/SMPrediction-dim19.txt"
smbin_corr_input = "validations/ATLAS/HIGG-2018-57/SMbin-corr.txt" 

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "validations/ATLAS/HIGG-2018-57/CgCGa-HIGG-2018-57-stxs.out"
outputplot = "validations/ATLAS/HIGG-2018-57/CgCGa-HIGG-2018-57-stxs.pdf"


# Scan ranges 
Cg_min = 0.8
Cg_max = 1.3
CGa_min = 0.8
CGa_max = 1.4

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CGa=1, Cg=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CV, CF"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="ff">1.</C>
  <C to="VV">1.</C>
  <C to="gammagamma">%(CGa)s</C>
  <C to="gg">%(Cg)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass': mass, 'CGa': CGa, 'Cg': Cg, 'precision': precision}
        
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

for Cg in np.linspace(Cg_min, Cg_max, grid_subdivisions):
    fresults.write('\n')
    for CGa in np.linspace(CGa_min, CGa_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, CGa=CGa, Cg=Cg, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            Cgmin = Cg
            CGamin = CGa
        fresults.write('%.5f    ' % Cg + '%.5f    ' % CGa + '%.5f     ' % m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at Cgluon, Cgamma, -2logL_min = ", Cgmin, CGamin, m2logLmin)

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

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
plt.tick_params(direction='in', labelsize=20, length=14, width=2)
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
Z = griddata((x, y), z2, (X,Y), method="linear")
# Import Official data from file 
dataload = open('validations/ATLAS/HIGG-2018-57/HIGG-2018-57-CgCga-Official.csv','r')
dorix = []
doriy = []
for line in dataload:
  fdat = line.split(',')
  dorix.append(float(fdat[0]))
  doriy.append(float(fdat[1]))
  
# Plotting the 68%, 95% and 99.7% CL contours
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

CS = ax.contour(xi,yi,Z,[2.3,5.99,11.83], colors=['silver'])
CS.levels = ['68% CL', '95% CL', '99.7% CL']
ax.clabel(CS, CS.levels, inline=1, fontsize=9, colors='k')

ax.set_aspect((0.33)/(CGa_max-CGa_min))

# best fit point
plt.plot([Cgmin],[CGamin], '*', c='w', ms=10)
# Standard Model 
plt.plot([1],[1], '+', c='k', ms=10)
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=3,c='b',marker='o',label='ATLAS official 68%')    
plt.legend(loc='lower right', scatterpoints = 3) 


# Title, labels, color bar...
plt.title("  Lilith-2.1, ATLAS HIGG-2020-16 validation" , fontsize=12, ha="center")
plt.xlabel(r'$C_g$',fontsize=25)
plt.ylabel(r'$C_\gamma$',fontsize=25)
plt.text(0.9, 1.3, r'Exp. input type = vn', fontsize=12, ha = 'left')
#plt.tight_layout()
fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/results")
print("***** done *****")

