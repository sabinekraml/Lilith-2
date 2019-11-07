###############################################################
#
# Lilith routine example (To run from /Lilith-1.1.x root folder)
#
# Constraints on the loop-induced reduced couplings CF, CV
# from a (CF, CV) fit
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
sys.path.append('.')
sys.path.append('../../../')
import lilith


# for plotting exp curves for comparison

######################################################################
# Parameters
######################################################################

print "***** reading parameters *****"

# Exprimental results
publication="ATLAS-HIGG-2018-57"
myexpinput ="my_exp.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Output
# Higgs mass to test
myhmass = 125.09
# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 50
# Output file
output = "output/CVCF_"+publication+".out"
# Output plot
outputplot = "plots/CVCF_"+publication+"_gamgam.pdf"
# Range of the scan
CF_min = 0.4
CF_max = 2.1
CV_min = 0.4
CV_max = 2.1

verbose = False
timer = False

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CF=1., CV=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings
       CF, CV"""

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="ff">%(CF)s</C>
  <C to="VV">%(CV)s</C>

  <extraBR>
    <BR to="invisible">0</BR>
    <BR to="undetected">0</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass': mass, 'CF': CF, 'CV': CV, 'precision': precision}

    return myInputTemplate % myInput

def getL(CF, CV):
    myXML_user_input = usrXMLinput(mass=myhmass, CF=CF, CV=CV, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l

######################################################################
# Scan initialization
######################################################################

print"***** scan initialization *****"

# Prepare output


# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False, timer=False)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

######################################################################
# Scan routine
######################################################################

m2logLmin = 10000
max = -1


print"***** running scan *****"

if not os.path.exists("output"):
	os.makedirs("output")

if os.path.isfile(output):
#	iscan='no'
	iscan='yes'
else:
	iscan='yes'

if iscan == "yes":
	fresults = open(output, 'w')
	for CF in np.linspace(CF_min, CF_max, grid_subdivisions):
		fresults.write('\n')
		for CV in np.linspace(CV_min, CV_max, grid_subdivisions):
			m2logL=getL(CF, CV)
			if m2logL < m2logLmin:
				m2logLmin = m2logL
				CFmin = CF
				CVmin = CV
			fresults.write('%.5f    ' % CV + '%.5f    ' % CF + '%.5f     ' % m2logL + '\n')

	fresults.close()

	print("m2logL min =", m2logLmin)
	print("CFmin, CVmin =", CFmin, CVmin)

	print"***** scan finalized *****"
	data = np.genfromtxt(output)
	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
else:
	fresults = open(output, 'r')
	data = np.genfromtxt(output)
	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	imin=np.where(z==z.min())[0][0]
	CVmin = x[imin]
	CFmin = y[imin]

######################################################################
# Plot routine
######################################################################


print "***** plotting *****"



# Preparing plot
fig = plt.figure(figsize=(8,7))
plt.subplots_adjust(left=0.115, bottom=0.11, right=0.97, top=0.90, wspace=0.4, hspace=0.3)
 
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)
plt.rc('axes', titlesize=22)     # fontsize of the axes title
plt.rc('axes', labelsize=17)
plt.rc('legend', fontsize=18)

ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

ax.tick_params(direction='in',length=6,width=2,which='major')
ax.tick_params(direction='in',length=4,width=1,which='minor')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.minorticks_on()



# Getting the data

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

ax.contourf(xi,yi,Z,[10**(-10),2.28,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

#ax.set_aspect((CF_max-CF_min)/(CV_max-CV_min))

#plt.plot([CVmin], [CFmin],'*', c='w', ms=10)
plt.plot([1],[1], '+', c='k', ms=10)

#  official ATLAS result
dt = np.dtype([('cx', float), ('cy', float)])
expCont = np.genfromtxt('data/fig12_CVCF_gammagamma.csv', dtype=dt)
plt.plot(expCont['cx'],expCont['cy'], '.', c='b', label='ATLAS Fig 12')
plt.legend(loc='lower right', ncol=1,frameon=False, prop={'size': 15},labelspacing=2.,handletextpad=0.5,fontsize=25)

# Title, labels, color bar...
plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion)+"     ", fontsize=14.5, ha="left")
plt.ylabel(r'$C_F$',fontsize=25)
plt.xlabel(r'$C_V$',fontsize=25)
plt.text(0.92, 1.35, "Data from "+publication, fontsize=13)


if not os.path.exists("plots"):
	os.makedirs("plots")
#plt.show()
plt.savefig(outputplot)

