###############################################################
#
# Lilith routine for (CV, CF) validation plots
#
# To run from /Lilith-1.x root folder
#
# Use the libraries matplotlib (plotting) and numpy (functions)
#
###############################################################

import sys, os
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from iminuit import Minuit

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print "***** reading parameters *****"

# Experimental results
exp_input = "data4validation/validation.list"
# Lilith precision mode
myprecision = "BEST-QCD"

# Higgs mass to test
myhmass = 125.09

# Output file
output = "validation/CMS/CgluCgam_2d.out"
# Output plot
outputplot = "validation/CMS/CgluCgam_2d.pdf"


# Range of the scan
Cg_min = 0.8
Cg_max = 1.6
CGa_min = 0.7
CGa_max = 1.4

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CGa=1, Cg=1, BRinv=0, BRund=0, precision="BEST-QCD"):
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
    <BR to="invisible">%(BRinv)s</BR>
    <BR to="undetected">%(BRund)s</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass': mass, 'CGa': CGa, 'Cg': Cg, 'BRinv': BRinv, 'BRund': BRund, 'precision': precision}
        
    return myInputTemplate%myInput

######################################################################
# * getL:        -2LogL for a given (CGa, Cg, BRinv) point
######################################################################

def getL(CGa, Cg, BRinv, BRund):
    myXML_user_input = usrXMLinput(mass=myhmass, CGa=CGa, Cg=Cg, BRinv=BRinv, BRund=BRund, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l


######################################################################
# Scan initialization
######################################################################

print "***** scan initialization *****"

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

print "***** running scan *****"

for CGa in np.linspace(CGa_min, CGa_max, grid_subdivisions):
    fresults.write('\n')
    for Cg in np.linspace(Cg_min, Cg_max, grid_subdivisions):
        #myXML_user_input = usrXMLinput(myhmass, CGa=CGa, Cg=Cg, precision=myprecision)
        #lilithcalc.computelikelihood(userinput=myXML_user_input)
        #m2logL = lilithcalc.l
        m = Minuit(getL, CGa=CGa, fix_CGa=True, Cg=Cg, fix_Cg=True, BRinv=0.1, limit_BRinv=(0,0.49), BRund=0.1, limit_BRund=(0,0.49), print_level=0, errordef=1, error_BRinv=0.05, error_BRund=0.05)
        m.migrad()
        m2logL = m.fval

        if m2logL < m2logLmin:
            m2logLmin = m2logL
            Cgmin = Cg
            CGamin = CGa
        fresults.write('%.5f    ' % CGa + '%.5f    ' % Cg + '%.5f     ' % m2logL + '\n')

fresults.close()

print "***** scan finalized *****"
print "minimum at Cgluon, Cgamma, -2logL_min = ", Cgmin, CGamin, m2logLmin


######################################################################
# Plot routine
######################################################################


print "***** plotting *****"

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=20, length=14, width=2)
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
Z = griddata(x, y, z2, xi, yi, interp="linear")

# Plotting the 68%, 95% and 99.7% CL contours

#ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
#              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

# Plotting the 1, 2 and 3 sigma contours

ax.contourf(xi,yi,Z,[10**(-10),2.3,6.18,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CGa_max-CGa_min)/(Cg_max-Cg_min))

plt.plot([CGamin],[Cgmin], '*', c='w', ms=10)
plt.plot([1],[1], '+', c='k', ms=10)

#  official CMS result
dt = np.dtype([('cx', float), ('cy', float)])
expCont = np.genfromtxt('validation/CMS/HIG-17-031_CgCGa-combined-2d-Grid.txt', dtype=dt)
plt.plot(expCont['cx'],expCont['cy'], '.', c='b', label='CMS official')
plt.legend(loc='lower right', fontsize=12)

# Title, labels, color bar...
plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion)+"     ", fontsize=14.5, ha="left")
plt.ylabel(r'$C_g$',fontsize=25)
plt.xlabel(r'$C_\gamma$',fontsize=25)
plt.text(0.75, 1.505, r'Data from CMS-HIG-17-031 & 17-023', fontsize=13)
plt.text(0.75, 1.455, r'BR($H\to$ inv.), BR($H\to$ undet.) profiled', fontsize=13)

#plt.tight_layout()
fig.set_tight_layout(True)
plt.show()
# Saving figure (.pdf)
#plt.savefig(outputplot)
fig.savefig("validation/CMS/HIG-17-031-CgCGa_BRinvBRund_profiled-v2.pdf")


