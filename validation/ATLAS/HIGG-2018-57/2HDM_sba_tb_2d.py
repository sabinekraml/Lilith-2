##################################################################
#
# Lilith routine example
#
# To put in Lilith-1.1.X/examples/python/
# To execute from /Lilith-1.X root folder
#
# Constraints on cos(beta-alpha) and tan(beta) in the 2HDM of
# Types I and II (in the sin(beta-alpha)>0 convention).
#
# Constraints are obtained through a likelihood-ratio test.
#
# Assumptions:
# * The lighter CP-even state h is identified with the observed
#   one (this fixes the coupling structure of the model)
#
# * The charged Higgs states are decoupled and thus
#   do not appear in the Higgs-photon-photon loop.
#   Only Higgs-fermion-fermion coupling modifications affect
#   this loop.
#
# Use the libraries matplotlib (plotting) and numpy (functions)
#
##################################################################

import sys, os
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('../../../')
import lilith

######################################################################
# Parameters
######################################################################

print "***** reading parameters *****"

# Experimental results

myexpinput ="my_exp.list"
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# 2HDM type = 1, 2
type = 2

# Number of steps for the square grid (cba,tb), need ~300 for fine grid but quite long to run
grid_subdivisions = 150

######################################################################

# Output file
if type == 1:
  output = "output/cba_tb_I_h_2d.out"
  outputplot = "plots/cba_tb_I_h_2d.pdf"
  cba_min = -0.6
  cba_max = 0.6
  tb_min = 0.1
  tb_max = 10.
  expfile='data/fig18a_2HDM_I_tanbeta_cosb-a-95CL.csv'
  textmodel="2HDM Type I"


if type == 2:
  output = "output/cba_tb_II_h_2d.out"
  outputplot = "plots/cba_tb_II_h_2d.pdf"
  cba_min = -0.2
  cba_max = 0.6
  tb_min = 0.1
  tb_max = 10.
  expfile='data/fig18b_2HDM_II_tanbeta_cosb-a-95CL.csv'
  textmodel="2HDM Type II"





######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, cba=0., tb=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings CGa, Cg"""
    
    sba = np.sqrt(1-cba**2)
    
    if type == 1:
      CV = sba
      CU = sba + cba/tb
      CD = sba + cba/tb
      CL = sba + cba/tb
    
    elif type == 2:
      CV = sba
      CU = sba + cba/tb
      CD = sba - cba*tb
      CL = sba - cba*tb

    
    elif type == 3:
      CV = sba
      CU = sba + cba/tb
      CD = sba + cba/tb
      CL = sba - cba*tb    

    elif type == 4:
      CV = sba
      CU = sba + cba/tb
      CD = sba - cba*tb
      CL = sba + cba/tb
    
    elif type == 5:
      MA= cba
      atb=np.sqrt( 1+ tb**2)
      afac=(MA**2 +MZ**2)*tb/(MZ**2 + MA**2*tb**2 - hmass**2*atb**2)
      su= 1/np.sqrt(1 +  afac**2)
      sd= su* afac
      CV = (sd+ su *tb)/atb
      CU= su*atb/tb
      CD= sd*atb
      CL= sd*atb
    else:
      print "Error: 2HDM type parameter should be 1 or 2"
      sys.exit()

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C>
  <C to="dd">%(CD)s</C>
  <C to="VV">%(CV)s</C>

  <C to="bb">%(CL)s</C> 
  <C to="tautau">%(CL)s</C>
  <C to="mumu">%(CL)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CV':CV, 'CU':CU, 'CD':CD, 'CL':CL, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

print "***** 2HDM scan initialization *****"

## Prepare output

#
## Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False,timer=False)
## Read experimental data
lilithcalc.readexpinput(myexpinput)
#
#
#######################################################################
## Scan routine
#######################################################################
#
m2logLmin=10000
max=-1

print "***** running 2HDM scan *****"

if not os.path.exists("output"):
	os.makedirs("output")

if os.path.isfile(output):
	iscan='no'
else:
	iscan='yes'

if iscan == "yes":
 fresults = open(output, 'w')
 for cba in np.linspace(cba_min, cba_max, grid_subdivisions):
    fresults.write('\n')
    for tb in np.linspace(tb_min, tb_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, tb=tb, cba=cba, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l                 #This is -2*Ln(L) at the (cba,tb) point
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            cbamin = cba
            tbmin = tb
        fresults.write('%.5f    '%cba +'%.5f    '%tb + '%.5f     '%m2logL + '\n')
 fresults.close()
 print "best fit, tbmin, cbamin, m2logLmin: ", tb, cba, m2logLmin
print "***** scan finalized *****"

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
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[el - z.min() for el in z]

print "2LogL Min: ", z.min()
# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z2, xi, yi, interp="linear")

levels = np.arange(-2.0, 1.601, 0.4)
cmap = matplotlib.cm.YlOrRd

# Plotting the 68%, 95% and 99.7% CL contours
ax.contourf(xi,yi,Z,[10**(-3),2.28,5.99,11.83],cmap=matplotlib.cm.get_cmap(cmap, len(levels) ))

# Title, labels, color bar etc.

dt = np.dtype([('cx', float), ('cy', float)])
expCont = np.genfromtxt(expfile, dtype=dt)
plt.plot(expCont['cx'],expCont['cy'], '.', c='b', label='ATLAS official 95% CL')


plt.title("  Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20.5, ha="left")

plt.ylabel(r'$\tan\beta$',fontsize=25)
plt.yscale('log')
if type == 2:
	plt.text(0.25, 0.32, textmodel, fontsize=13)
	plt.text(0.25, 0.25, "Data from ATLAS-HIGG-2018-57", fontsize=13)

if type == 1:
	plt.text(0.15, 0.32, textmodel, fontsize=13)
	plt.text(0.15, 0.25, "Data from", fontsize=13)
	plt.text(0.15, 0.20, "ATLAS-HIGG-2018-57", fontsize=13)

plt.legend(bbox_to_anchor=(0.97,0.01),loc='lower right', ncol=1,frameon=True, prop={'size': 13},labelspacing=0.2,handletextpad=0.5,fontsize=15)
if type == 1:
  plt.xlim([cba_min-0.01,cba_max +0.01])
  plt.xlabel(r'$\cos(\beta-\alpha)$',fontsize=25)
if type == 2:
  plt.xlim([cba_min-0.01,cba_max +0.01])
  plt.xlabel(r'$\cos(\beta-\alpha)$',fontsize=25)

plt.tight_layout()

if not os.path.exists("plots"):
	os.makedirs("plots")
#plt.show()
plt.savefig(outputplot)

