##################################################################
#
# Lilith routine example
#
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
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
# Uses the libraries matplotlib (plotting) and numpy (functions)
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from math import floor, log, sqrt, sin, cos, atan
from cmath import sqrt as csqrt
from cmath import asin as casin

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/cba_tb/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

#resultstype = "latestRun2.list"
resultstype = "thisRun2.list"

# 2HDM type = 1, 2
yukawatype = "I"

# Experimental results
exp_input = lilith_dir+"validations/STU/finalanalysis/" + resultstype

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

mHpm = 500

# Output files
if resultstype == "latestRun2.list":
	output = validation_dir+"cba_tb_" + yukawatype + "_2d_36fb_Hpm.out"
	outputplot = validation_dir+"cba_tb_" + yukawatype + "_2d_36fb_Hpm.pdf"
if resultstype == "thisRun2.list":
	output = validation_dir+"cba_tb_"+ yukawatype + "_2d_140fb_Hpm.out"
	outputplot = validation_dir+"cba_tb_" + yukawatype + "_2d_140fb_Hpm.pdf"


# Scan ranges
if yukawatype == "I":
  cba_min = -0.5
  cba_max = 0.5
  tb_min = 0.1
  tb_max = 10.

if yukawatype == "II":
  cba_min = -0.1
  cba_max = 0.6
  tb_min = 0.1
  tb_max = 10.

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 150

########################################################################
# SM and H+ LO contribution to the reduced H-gamma-gamma coupling CGa
########################################################################

# Values

mW = 80.398 
v = 246
mZ = 91.1876

mt = 173.1
mb = 4.75
mc = 1.4
mtau = 1.777
sW2 = 0.23116

# New values ????

def fhiggs(t):
    if t<=1.:
        return casin(sqrt(t))**2.
    else:
        return -(log((sqrt(t)+sqrt(t-1.))/(sqrt(t)-sqrt(t-1.)))-np.pi*1j)**2./4.

def A0(tau):
    return -1./tau *(1.-1./tau * fhiggs(tau))

def A12(tau):
    return 2./tau *(1.+(1.-1./tau) * fhiggs(tau))

def A1(tau):
    return -(3.*tau+2.*tau**2. +3.*(2.*tau-1.) * fhiggs(tau))/tau**2


#def get_CGa(a, tb, b, sinba, cosba, l1, l2, l3, l4, l5):
def get_CGa(a, tb, b, sinba, cosba):
	""" 
      Returns CGa computed from the SM particles contribution plus 
      Hpm contribution.
	"""

#	if yukawatype == "I":
#		CV = sinba
#		CU = np.cos(a)/np.sin(b)
#		CD = np.cos(a)/np.sin(b)
#    
#	elif yukawatype == "II":
#		CV = sinba
#		CU = np.cos(a)/np.sin(b)
#		CD = -np.sin(a)/np.cos(b)

	if yukawatype == "I":
		CV = sinba
		CU = sinba + cosba/tb
		CD = sinba + cosba/tb
    
	elif yukawatype == "II":
		CV = sinba
		CU = sinba + cosba/tb
		CD = sinba - cosba*tb

	A12t = A12((hmass/(2.*mt))**2.)
	A12c = A12((hmass/(2.*mc))**2.)
	A12b = A12((hmass/(2.*mb))**2.)
	A12tau = A12((hmass/(2.*mtau))**2.)
	A1W = A1((hmass/(2.*mW))**2.)
	A0Hpm = A0((hmass/(2.*mHpm))**2.)

#	Z3 = 0.25*np.sin(2*b)**2 * (l1 + l2 + -2*( l3 + l4 + l5 ) ) + l3
#	Z7 = -0.5*np.sin(2*b) * ( l1*np.sin(b)**2 - l2*np.cos(b)**2 + ( l3 + l4 + l5 )*np.cos(2*b) ) 
#	ChHpmHpm = -v*(Z3*sinba + Z7*cosba)

	m122 = np.cos(a)**2*hmass*2/tb
	mbarre2 = m122/(np.sin(b)*np.cos(b))
	ChHpmHpm = (-1/v)*( (hmass**2 + 2*(mHpm**2 - mbarre2) )*sinba + 2*(1/np.tan(2*b))*(hmass**2-mbarre2)*cosba )

	SM_particles_amplitude = CU*4./3.*(A12t + A12c) + CD*1./3.*A12b + CD*A12tau + CV*A1W
	SM_amplitude = 4./3.*(A12t + A12c) + 1./3.*A12b + A12tau + A1W
	Hpm_amplitude = ChHpmHpm*A0Hpm
  
	CGa = sqrt( abs(SM_particles_amplitude + Hpm_amplitude)**2./abs(SM_amplitude)**2. )

#	CGa = np.sqrt( abs ((CU*4./3.*(A12t + A12c) + CD*1./3.*A12b + CD*A12tau + CV*A1W ) / ( 4./3.*(A12t + A12c) + 1./3.*A12b + A12tau + A1W ) ))

	C = [CV, CU, CD, CGa]

	return C

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CV=1, CU=1, CD=1, CGa=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CU, CD, CV, CGa"""

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C>
  <C to="dd">%(CD)s</C>
  <C to="VV">%(CV)s</C>
  <C to="gammagamma">%(CGa)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CV':CV, 'CU':CU, 'CD':CD, 'CGa':CGa, 'precision':precision}
        
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

print("***** running 2HDM scan *****")

for cba in np.linspace(cba_min, cba_max, grid_subdivisions):
    fresults.write('\n')
    for tb in np.linspace(tb_min, tb_max, grid_subdivisions):
        b = np.arctan(tb)
        a = b - np.arccos(cba)
        sinba = np.sqrt(1-cba**2) 	
        C = get_CGa(a=a, tb=tb, b=b, sinba=sinba, cosba=cba)
        myXML_user_input = usrXMLinput(mass=hmass, CV=C[0], CU=C[1], CD=C[2], precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l                 #This is -2*Ln(L) at the (cba,tb) point
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            cbamin = cba
            tbmin = tb
        fresults.write('%.5f    '%cba +'%.5f    '%tb + '%.5f     '%m2logL + '\n')
fresults.close()

print("***** scan finalized *****")

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure(figsize=(7.5,7))
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
z2=[el - z.min() for el in z]

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)

try:
# griddata using Natural Neighbor (nn) interpolation
    Z = griddata((x, y), z2, (X, Y), method="nearest")

except:
# If you don't have natgrid installed, use instead linear interpolation
    Z = griddata((x, y), z2, (X, Y), method="linear")

levels = np.arange(-2.0, 1.601, 0.4)
cmap = matplotlib.cm.YlOrRd

# Plotting the 68%, 95% and 99.7% CL contours (filled regions)
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],cmap=matplotlib.cm.get_cmap(cmap, len(levels) ))

#ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99,11.83],colors=['#ff3300','#ffa500','#ffff00'], \
#              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

# Title, labels, color bar etc.
plt.title("                              Lilith-2.1", fontsize=15, ha="left")
plt.xlabel(r'$\cos(\beta-\alpha)$',fontsize=25)
plt.ylabel(r'$\tan\beta$',fontsize=25)
plt.yscale('log')
if resultstype == "latestRun2.list":
  plt.xlim([-0.601,0.601])
  plt.text(-0.55, 0.2, '2HDM Type-' + yukawatype, fontsize=15)
  plt.text(-0.55, 0.16, r'(Run 2, 36$fb^{-1}$, ATLAS+CMS)', fontsize=12)
if resultstype == "thisRun2.list":
  plt.xlim([-0.601,0.601])
  plt.text(-0.55, 0.2, '2HDM Type-' + yukawatype, fontsize=15)
  plt.text(-0.55, 0.16, r'(Run 2, ATLAS 36$fb^{-1}$ + CMS 140$fb^{-1}$)', fontsize=12)


fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", validation_dir)
print("***** done *****")
