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

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
#lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "validations/ATLAS/HIGG-2018-28/latestRun2.list"
# Sm predictions     
smpred_input = "validations/ATLAS/HIGG-2018-28/SMbin-prediction-noThunc.txt"
smbin_corr_input = "validations/ATLAS/HIGG-2018-28/SMbin-corr.txt"
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/CVCF_2d.out"
#outputplot = "validations/ATLAS/HIGG-2018-28/CVCF_2d_Fig10a_no-theo-unc.pdf"
outputplot = "validations/ATLAS/HIGG-2018-28/CVCF_2d_STXS.pdf"

# Scan ranges
CV_min = 0.85
CV_max = 1.2
CF_min = 0.4
CF_max = 2.6

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100


######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125, CV=1, CF=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CV, CF"""

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="tt">%(CF)s</C>
  <C to="bb">%(CF)s</C>
  <C to="cc">%(CF)s</C>
  <C to="tautau">%(CF)s</C>
  <C to="ZZ">%(CV)s</C>
  <C to="WW">%(CV)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass': mass, 'CV': CV, 'CF': CF, 'precision': precision}

    return myInputTemplate % myInput


######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

# Prepare output
fresults = open(output, 'w')

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False, timer=False)
# Read experimental data
lilithcalc.readexpinput(exp_input)

# Read SM prediction input and correlation 

lilithcalc.readsmpred(smpred_input)

lilithcalc.readsmcorr(smbin_corr_input)
######################################################################
# Scan routine
######################################################################

m2logLmin = 10000
max = -1

print("***** running scan *****")

for CV in np.linspace(CV_min, CV_max, grid_subdivisions):
    fresults.write('\n')
    for CF in np.linspace(CF_min, CF_max, grid_subdivisions):
#        print("testing...") 
        myXML_user_input = usrXMLinput(hmass, CV=CV, CF=CF, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CVmin = CV
            CFmin = CF
        fresults.write('%.5f    ' % CV + '%.5f    ' % CF + '%.5f     ' % m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CV, CF, -2logL_min = ", CVmin, CFmin, m2logLmin)

######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 10
matplotlib.rcParams['ytick.major.pad'] = 10
#plt.locator_params(axis='y', nbins=10)
#plt.locator_params(axis='x', nbins=10)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=10, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(output)

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2 = []
for z_el in z:
    z2.append(z_el - z.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")

# Plotting the 68%, 95% and 99.7% CL regions
ax.contourf(xi, yi, Z, [10 ** (-10), 2.3, 5.99, 11.83], colors=['#ff3300', '#ffa500', '#ffff00'], \
            vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((CV_max - CV_min) / (CF_max - CF_min))

# best fit point
plt.plot([CVmin], [CFmin], '*', c='w', ms=10)

# Standard Model 
plt.plot([1], [1], '+', c='k', ms=10, label="SM")

# Official best fit point
plt.plot([1.022857143], [0.878732748], '*', c='b', ms=10, label="Official best fit")

# read data for official 68% and 95% CL contours & plot (added by TQL)
expdata = np.genfromtxt('validations/ATLAS/HIGG-2018-28/HIGG-2018-28-CVCF-Grids.txt')
xExp = expdata[:, 0]
yExp = expdata[:, 1]
plt.plot(xExp, yExp, '.', markersize=4, color='blue', label="ATLAS official")
#plt.scatter(dorix,doriy,s=6,c='b',marker='o',label='ATLAS official')    
plt.legend(loc='best', scatterpoints = 3) 

# Title, labels, color bar...
plt.title("  Lilith - STXS - test - HIGG-2018-28", fontsize=14, ha="center")
plt.xlabel(r'$C_V$', fontsize=12)
plt.ylabel(r'$C_F$', fontsize=12)
#plt.text(0.6, 2.15, r'Exp. input: ATLAS-HIGG-2018-28_Fig10a_no-theo-unc.xml', fontsize=12)
#plt.text(0.87, 2.2, r'ggF theo. correlation', fontsize=10)
#plt.text(0.7, 2.0, r'very fisrt test, no theoretical uncertainty', fontsize=8)
#plt.text(0.7, 2.0, r'very fisrt test, with theo. uncertainty, approx. 3', fontsize=8)
plt.text(0.87, 2.3, r'no theo. uncertainty', fontsize=10)
plt.text(0.87, 2.4, r'type = vn', fontsize=10)




fig.set_tight_layout(True)

# plt.show()
#set aspect ratio to 0.8
ratio = 0.8
x_left, x_right = ax.get_xlim()
y_low, y_high = ax.get_ylim()
ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)


# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/validations/ATLAS/HIGG-2018-28/")
print("***** done *****")

