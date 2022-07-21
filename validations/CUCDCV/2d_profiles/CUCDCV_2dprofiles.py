##################################################################
#
# Lilith routine for (CU, CD) validation plots
#
##################################################################
import time
import sys, os
from scipy.interpolate import griddata
from iminuit.minimize import minimize
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

if (not os.path.exists("results")):
    os.mkdir("results")

if (not os.path.exists("plots")):
    os.mkdir("plots")
    
validation_dir = lilith_dir+"validations/CUCDCV/2d_profiles/"

######################################################################
#Execution time
start_time = time.time()
######################################################################

#######################################################################
# Parameters
######################################################################

# Exprimental results
myexpinput = lilith_dir+"data/latestRun2_140fb-1.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Output
output = validation_dir+"results/CUCD-CV_2dprofiles_140fb-1_minuit.out"
outputplot = "plots/CUCD-CV_2dprofiles_140fb-1_minuit.pdf"

# Scan ranges
CU_min = 0.6
CU_max = 1.4
CD_min = 0.6
CD_max = 1.4
CV_min = 0.5
CV_max = 1.5

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

# Higgs mass to test
hmass = 125.09

verbose=False
timer=False

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose, timer)

# Read experimental data
lilithcalc.readexpinput(myexpinput)

print("\nTask: profile-likelihood analysis in (CU,CD,CV) model")
print("Experimental input:", myexpinput)

######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given (CU, CD, CV) point
######################################################################

def usrXMLinput(mass=125.09, CU=1, CD=1, CV=1, precision="BEST-QCD"):
	"""generate XML input from reduced couplings CU, CD"""
    
	myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C> 
  <C to="dd">%(CD)s</C>
  <C to="WW">%(CV)s</C>
  <C to="ZZ">%(CV)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
	myInput = {'mass':mass, 'CU':CU, 'CD':CD, 'CV':CV, 'precision':precision}
        
	return myInputTemplate%myInput


######################################################################
# * getL:        -2LogL for a given (CU, CD, CV) point
######################################################################

def getL(CV, CU, CD):
    myXML_user_input = usrXMLinput(mass=hmass, CU=CU, CD=CD, CV=CV[0], precision=my_precision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l


######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

# Start values for fit parameters
CU0 = 1
CD0 = 1
CV0 = 1


# Prepare output
fresults = open(output, 'w')

######################################################################
# Scan routine
######################################################################

m2logLmin=10000
max=-1

print("***** running scan *****")

for CU in np.linspace(CU_min, CU_max, grid_subdivisions):
    fresults.write('\n')
    for CD in np.linspace(CD_min, CD_max, grid_subdivisions):
        getLminimized = minimize(getL, [1] , args=(CU, CD), method='migrad', bounds = [0.5,1.5])
        m2logL = getLminimized.fun
        CV = getLminimized.x[0]
        
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            CUmin = CU
            CDmin = CD
        fresults.write('%.5f    '%CU +'%.5f    '%CD + '%.5f     '%m2logL +'%.5f     '%CV + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at CU, CD, -2logL_min = ", CUmin, CDmin, m2logLmin)

print("--- %s seconds ---" % (time.time() - start_time))
