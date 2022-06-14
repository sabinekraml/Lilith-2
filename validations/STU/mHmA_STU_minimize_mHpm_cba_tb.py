##################################################################
#
# (S,T,U) using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 2d likelihood contour on the $m_H$, $m_A$ plane with $m_Hpm$, $\cos(\beta-\alpha)$, $\tan(\beta)$ minimized at each point
#
##################################################################

import sys, os
from scipy.interpolate import griddata
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/rangeminimize_cba_tb/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Values
Scen = 0.06
Ssigma = 0.10
Tcen = 0.11
Tsigma = 0.12
Ucen = 0.14
Usigma = 0.09
STcorrelation = 0.9
SUcorrelation = -0.59
TUcorrelation = -0.85
CEN_STU = np.array([Scen, Tcen, Ucen])
SIG_STU = np.diag([Ssigma, Tsigma, Usigma])
COR_STU = np.array(([1, STcorrelation, SUcorrelation],
                [STcorrelation, 1 , TUcorrelation],
                [SUcorrelation, TUcorrelation, 1]))	
C_STU = SIG_STU.dot(COR_STU).dot(SIG_STU)
C_STU_inv = np.linalg.inv(C_STU)

# Scan ranges
mA_min = 200
mA_max = 2000
mH_min = 200
mH_max = 2000
mHpm_min = 200
mHpm_max = 2000

# Experimental results
#exp_input = lilith_dir+"validations/STU/" + "thisRun2.list"
exp_input = lilith_dir+"validations/STU/" + "latestRun2.list"


# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
my_hmass = 125

# 2HDM type = 1, 2
type = 1

# Scan ranges
if type == 1:
  cba_min = -0.25
  cba_max = 0.25
  tb_min = 0.1
  tb_max = 10

if type == 2:
  cba_min = -0.05
  cba_max = 0.05
  tb_min = 0.1
  tb_max = 10

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 50

# Output files
if type == 1:
  output = validation_dir+"mHmA_STU_mHpm_cba_tb_" + str(grid_subdivisions) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + "_" + str(cba_min) + "-" + str(cba_max) + "_" + str(tb_min) + "-" + str(tb_max) + "_I" + ".out"
  outputplot = validation_dir+"mHmA_STU_mHpm_cba_tb" + str(grid_subdivisions) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + "_" + str(cba_min) + "-" + str(cba_max) + "_" + str(tb_min) + "-" + str(tb_max) + "_I" + ".pdf"

if type == 2:
  output = validation_dir+"mHmA_STU_mHpm_cba_tb_" + str(grid_subdivisions) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + "_" + str(cba_min) + "-" + str(cba_max) + "_" + str(tb_min) + "-" + str(tb_max) + "_II" + ".out"
  outputplot = validation_dir+"mHmA_STU_mHpm_cba_tb" + str(grid_subdivisions) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + "_" + str(cba_min) + "-" + str(cba_max) + "_" + str(tb_min) + "-" + str(tb_max) + "_II" + ".pdf"

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

def func(X, mH, mA):
		mHpm, cba, tb = X[0], X[1], X[2]
		sinba = np.sqrt(1-cba**2)
		S, T, U = calc.Scalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = sinba), calc.Tcalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = sinba), calc.Ucalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = sinba)
		X_STU = [S, T, U]
#		print("X = ", X)
		L2t_STU = C_STU_inv.dot(X_STU-CEN_STU).dot((X_STU-CEN_STU).T)

		myXML_user_input = usrXMLinput(mass=my_hmass, cba=cba, tb=tb, precision=my_precision)
		lilithcalc.computelikelihood(userinput=myXML_user_input)
		L2t_cba_tb = lilithcalc.l

		L2t = L2t_STU + L2t_cba_tb
#		print("blablabla = ", mHpm, cba, tb, mH, mA, L2t)

		return L2t

######################################################################
# Scan initialization
######################################################################

m2logLmin=10000
cba0=0.001
tb0=1.001
i=1

print("***** running scan *****")

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    if i==1:
        print("mH = ", mH)
    if i%10==0:
        print("mH = ", mH)
    i+=1
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        funcminimized = minimize(func, [(mH+mA)/2,cba0,tb0] , args=(mH, mA), method='SLSQP', bounds=((mHpm_min,mHpm_max),(cba_min,cba_max),(tb_min,tb_max)), options={'ftol': 1e-3, 'eps': 0.1} )
        m2logL = funcminimized.fun
        fit = funcminimized.x
        if funcminimized.success == False :
            print("Could not minimize for (mH, mA) = ", mH, mA)
#        print("mH, mA, mHpm, m2logLminpas = ", mH, mA, mHpmfit, m2logL)
#        m2logLminpas=10000
#        for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
#            m2logLpas = func(mH=mH, mA=mA, mHpm=mHpm)
#            if m2logLpas < m2logLminpas:
#                m2logLminpas = m2logLpas
#                mHpmminpas = mHpm
#        print("mH, mA, mHpm, m2logLminpas = ", mH, mA, mHpmminpas, m2logLminpas, "\n")
        fresults.write('%.5f    '%mH + '%.5f    '%mA + '%.5f    '%m2logL + '%.5f    '%fit[0] + '%.5f    '%fit[1] + '%.5f    '%fit[2] + '\n')
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            mHmin = mH
            mAmin = mA
            mHpmmin = fit[0]
            cbamin = fit[1]
            tbmin = fit[2]

fresults.write('\n' + '%.5f    '%mHmin + '%.5f    '%mAmin + '%.5f    '%m2logLmin + '%.5f    '%mHpmmin + '%.5f    '%cbamin + '%.5f    '%tbmin)
fresults.close()

print("***** scan finalized *****")
print("minimum at mH, mA, mHpm, cba, tb -2logL_min = ", mHmin, mAmin, mHpmmin, cbamin, tbmin, m2logLmin)

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

x = data[0:-1,0]
y = data[0:-1,1]
z = data[0:-1,2]


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
ax.contourf(xi,yi,Z,[10**(-10),3.506,7.815],colors=['#ff3300','#ffa500',])

ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

# best fit point
plt.plot([data[-1,0]],[data[-1,1]], '+', markersize=8, color = 'black', label = 'best fit')
plt.legend(loc='upper right')

# Title, labels, color bar...
#plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$m_H$[GeV]',fontsize=16)
plt.ylabel(r'$m_A$[GeV]',fontsize=16)
#plt.text(mH_min + 100, mA_max - 150, r'Contour plot in the $m_H$, $m_A$ plane with $m_{H^{\pm}}$ minimized at each point', fontsize=7)
#plt.text(mH_min + 100, mA_max - 250, fr"Range of $m_{{H^{{\pm}}}}$ = ({mHpm_min},{mHpm_max}), with $\cos(\beta - \alpha) = 0$", fontsize=7)
#plt.text(mH_min + 100, mA_max - 350, fr"Best point ($m_H, m_A$) = ({data[-1,0]:.0f}, {data[-1,1]:.0f}) with $\chi^{2}$ = {data[-1,2]:.3f} and $m_{{H^{{\pm}}}}$ = {data[-1,3]:.0f}", fontsize=7) 
#plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}", fontsize=9)
#plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
#plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", validation_dir)
print("***** done *****")
