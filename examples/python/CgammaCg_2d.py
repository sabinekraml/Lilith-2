###############################################################
#
# Lilith routine example (To run from /Lilith-1.1.x root folder)
#
# Constraints on the loop-induced reduced couplings Cgamma, Cg
# from a (Cgamma, Cg) fit
#
# Use the libraries matplotlib (plotting) and numpy (functions)
#
###############################################################

import sys, os
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print "***** reading parameters *****"

# Experimental results
exp_input = "data/latest.list"
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.

# Output file
output = "results/CgammaCg_2d.out"
# Output plot
outputplot = "results/CgammaCg_2d.pdf"


# Range of the scan
CGa_min = 0.
CGa_max = 2.
Cg_min = 0.
Cg_max = 2.

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125., CGa=1, Cg=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings
       CGa, Cg"""
    
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
    myInput = {'mass':mass, 'CGa':CGa, 'Cg':Cg, 'precision':precision}
        
    return myInputTemplate%myInput

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
        myXML_user_input = usrXMLinput(hmass, CGa=CGa, Cg=Cg, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CGamin = CGa
            Cgmin = Cg
        fresults.write('%.5f    '%CGa +'%.5f    '%Cg + '%.5f     '%m2logL + '\n')

fresults.close()

print "***** scan finalized *****"

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
ax.contour(xi,yi,Z,[2.3,5.99,11.83],linewidths=[2.5,2.5,2.5],colors=["#B22222","#FF8C00","#FFD700"])
cax=ax.imshow(Z, vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], \
              aspect=(CGa_max-CGa_min)/(Cg_max-Cg_min), cmap=plt.get_cmap("rainbow"))

# Title, labels, color bar...
plt.title("  Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=14.5, ha="left")
plt.xlabel(r'$C_\gamma$',fontsize=25)
plt.ylabel(r'$C_g$',fontsize=25)

plt.tight_layout()
cbar = fig.colorbar(cax)
cbar.set_label(r"$\Delta(-2\log L)$",fontsize=25)

# Saving figure (.pdf)
plt.savefig(outputplot)


