###############################################################
#
# Lilith routine for (C_gluon, C_gamma) validation plots
#
###############################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = "/Users/kraml/Lilith-dev/william/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CMS/HIG-19-015-sk/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = validation_dir+"thisRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
output = validation_dir+"HIG-19-015-CgCga_2d.out"
outputplot = validation_dir+"HIG-19-015-CgCga_2d.pdf"

# Scan ranges 
Cg_min = 0.6
Cg_max = 1.4
CGa_min = 0.7
CGa_max = 1.5

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

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
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

# Plotting the 68%, 95% and 99.7% CL contours
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

# to add contour lines and labels
#CS = ax.contour(xi,yi,Z,[2.3,5.99,11.83], colors=['silver'])
#CS.levels = ['68% CL', '95% CL', '99.7% CL']
#ax.clabel(CS, CS.levels, inline=1, fontsize=9, colors='k')

ax.set_aspect((0.6)/(0.6))

# best fit point
plt.plot([Cgmin],[CGamin], '*', c='w', ms=10)
# Standard Model 
plt.plot([1],[1], '+', c='k', ms=10)

# read data for official 68% and 95% CL contours & plot  
expdata = np.genfromtxt('HIG-19-015_CgCGa-Grid.txt')
xExp = expdata[:,0]
yExp = expdata[:,1]
plt.plot(xExp,yExp,'.',markersize=4, color = 'blue', label="CMS official")
plt.legend()

# Title, labels, color bar...
plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$C_g$',fontsize=18)
plt.ylabel(r'$C_\gamma$',fontsize=18)
plt.text(0.68, 0.8, r'Data from CMS-HIG-19-015', fontsize=12)
plt.text(0.68, 0.75, r'(Fig. 16 + Aux. Fig. 3)', fontsize=11)

#plt.tight_layout()
fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", validation_dir)

