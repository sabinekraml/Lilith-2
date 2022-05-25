##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 2d likelihood contour on the $m_H$, $m_A$ plane with $m_Hpm$ minimized at each point
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

#validation_dir = lilith_dir+"validations/STU/"
validation_dir = lilith_dir+"validations/STU/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Values
Scen = 0.15
Ssigma = 0.08
Tcen = 0.27
Tsigma = 0.06
STcorrelation = 0.93

# Scan ranges
#mA_min = 600
#mA_max = 1400
#mH_min = 600
#mH_max = 1400
#mHpm_min = 500
#mHpm_max = 1000
mA_min = 200
mA_max = 400
mH_min = 200
mH_max = 400
mHpm_min = 200
mHpm_max = 400

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 10

# Experimental results
exp_input = validation_dir+"thisRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
my_hmass = 125

# 2HDM type = 1, 2
type = 1

# Output files
if type == 1:
  output = validation_dir+"mHmA_ST_mHpm_cba_tb_I.out"
  outputplot = validation_dir+"mHmA_ST_mHpm_cba_tb_I.pdf"

if type == 2:
  output = validation_dir+"mHmA_ST_mHpm_cba_tb_II.out"
  outputplot = validation_dir+"mHmA_ST_mHpm_cba_tb_II.pdf"

# Scan ranges
if type == 1:
  cba_min = -0.25
  cba_max = 0.25
  tb_min = 0.1
  tb_max = 10.

if type == 2:
  cba_min = -0.05
  cba_max = 0.05
  tb_min = 0.1
  tb_max = 10.

#Precision
#mHpm_precision = 30
#cba_precision = 15
#tb_precision = 15
mHpm_precision = 10
cba_precision = 10
tb_precision = 10

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
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125, cba=0., tb=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings CU, CD, CV"""
    
    sba = np.sqrt(1-cba**2)
    
    if type == 1:
      CV = sba
      CU = sba + cba/tb
      CD = sba + cba/tb
    
    elif type == 2:
      CV = sba
      CU = sba + cba/tb
      CD = sba - cba*tb

    else:
      print("Error: 2HDM type parameter should be 1 or 2")
      sys.exit()

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C>
  <C to="dd">%(CD)s</C>
  <C to="VV">%(CV)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CV':CV, 'CU':CU, 'CD':CD, 'precision':precision}
        
    return myInputTemplate%myInput


######################################################################
# Likelihood Calculation
######################################################################

def func(mH, mA, mHpm, cba):
		sinba = np.sqrt(1-cba**2)
		z10, z20 = Scen, Tcen
		sig1m, sig1p = Ssigma, Ssigma
		sig2m, sig2p = Tsigma, Tsigma
		p = STcorrelation

		z1, z2 = calc.Scalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1), calc.Tcalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1)

		V1 = sig1p * sig1m
		V1e = sig1p - sig1m
		V2 = sig2p * sig2m
		V2e = sig2p - sig2m
		V1f = V1 + V1e * (z1 - z10)
		V2f = V2 + V2e * (z2 - z20)
		L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#		print("mH, mA, mHpm, L2t = ", mH, mA, mHpm, L2t)
		return L2t

def getL(cba, tb):
    myXML_user_input = usrXMLinput(mass=my_hmass, cba=cba, tb=tb, precision=my_precision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l

######################################################################
# Scan initialization
######################################################################

m2logLmin=10000

print("***** running scan *****")

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logLmincurrent=10000
        for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
            for cba in np.linspace(cba_min, cba_max, cba_precision):
                for tb in np.linspace(tb_min, tb_max, tb_precision):
                    m2logL = func(mH=mH, mA=mA, mHpm=mHpm, cba=cba) + getL(cba=cba, tb=tb)
                    if m2logL < m2logLmincurrent:
                        m2logLmincurrent = m2logL
                    print("mH, mA, mHpm, cba, tb, m2logL, m2logLmincurrent = ", mH, mA, mHpm, cba, tb, m2logL, m2logLmincurrent)
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logLmincurrent + '\n')
        if m2logLmincurrent < m2logLmin:
            m2logLmin = m2logLmincurrent
            mHmin = mH
            mAmin = mA

fresults.close()

print("***** scan finalized *****")
print("minimum at mH, mA, -2logL_min = ", mHmin, mAmin, m2logLmin)

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

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
Z = griddata((x, y), z2, (X, Y), method="linear")


# Plotting the 68% and 95% CL regions
#ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99, 7.81],colors=['#ff3300','#ffa500', '#ffff00'])#, \
#              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
ax.contourf(xi,yi,Z,[10**(-10),3.506,7.815],colors=['#ff3300','#ffa500',])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

#sc = ax.scatter(x, y, c=z2)
#cbar = fig.colorbar(sc)

ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

# best fit point
plt.plot([mHmin],[mAmin], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$m_H$[GeV]',fontsize=16)
plt.ylabel(r'$m_A$[GeV]',fontsize=16)
plt.text(mH_min + 100, mA_max - 150, r'Contour plot in the $m_H$, $m_A$ plane with $m_{H^{\pm}}$ minimized at each point', fontsize=7)
plt.text(mH_min + 100, mA_max - 250, fr"Range of $m_{{H^{{\pm}}}}$ = ({mHpm_min},{mHpm_max}), with $\cos(\beta - \alpha) = 0$", fontsize=7)
plt.text(mH_min + 100, mA_max - 350, f"Best point ($m_H, m_A$) = ({mHmin:.0f}, {mAmin:.0f})", fontsize=7)
#plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}", fontsize=9)
#plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
#plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "validations/STU/")
print("***** done *****")
