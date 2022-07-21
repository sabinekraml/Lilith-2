##################################################################
#
# Lilith routine for (CU, CV) validation plots
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CUCDCV/2d_simple/"

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = lilith_dir + "data/latestRun2_140fb-1.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
if (not os.path.exists("plots")):
    os.mkdir("plots")

if (not os.path.exists("results")):
    os.mkdir("results")
    
output = validation_dir+"results/CVCD_140fb-1.out"
outputplot = validation_dir+"plots/CVCD_140fb-1.pdf"


# Scan ranges
CD_min = 0.76
CD_max = 1.28
CV_min = 0.9
CV_max = 1.22

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100


######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CV=1, CD=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings CV, CU"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>


  <C to="WW">%(CV)s</C>
  <C to="ZZ">%(CV)s</C>
  <C to="dd">%(CD)s</C> 
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CV':CV, 'CD':CD, 'precision':precision}
        
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

for CV in np.linspace(CV_min, CV_max, grid_subdivisions):
    fresults.write('\n')
    for CD in np.linspace(CD_min, CD_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, CV=CV, CD=CD, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CVmin = CV
            CDmin = CD
        fresults.write('%.5f    '%CV +'%.5f    '%CD + '%.5f     '%m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CV, CD, -2logL_min = ", CVmin, CDmin, m2logLmin)


