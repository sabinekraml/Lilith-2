###############################################################
#
# Lilith routine example
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Constraints on reduced couplings in the (CV, CF) model. 
# CV: common coupling to EW gauge bosons,
# CF: common coupling for fermions.
# The loop-induced effective couplings CGa and Cg are
# determined from CV and CF.
#
# 1-dimensional likelihood profiles are obtained from a
# profile-likelihood analysis
#
# Use the libraries matplotlib, iminuit
#
###############################################################

import sys, os
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CUCDCV/1d_profiles"

if (not os.path.exists("results")):
    os.mkdir("results")
    
import iminuit

######################################################################
# Parameters
######################################################################

# Exprimental results
myexpinput = lilith_dir+"data/latestRun2.list"

# Lilith precision mode
myprecision = "BEST-QCD"


# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

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
    
def getL(CU, CD, CV):
    myXML_user_input = usrXMLinput(mass=myhmass, CU=CU, CD=CD, CV=CV, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l
    
    
######################################################################
# Calculations
######################################################################

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose, timer)

# Read experimental data
lilithcalc.readexpinput(myexpinput)

print("\n***** performing model fit: migrad, hesse, minos *****")

# Start values for fit parameters
CU0 = 1
CD0 = 1
CV0 = 1

# Initialize the fit; parameter starting values and limits

m = Minuit(getL, CU=CU0, CD=CD0, CV=CV0)

m.limits = [(0,3), (0,3),(0,3)]

# Minimization and error estimation
m.migrad()
m.hesse()   # run covariance estimator
m.minos()


print("\nbest-fit point:", m.values) 
print("\nHesse errors:", m.errors)
print("\nMinos errors:")
for key, value in list(m.merrors.items()):
    print(key, value)



print("\n***** getting the 1-dimensional likelihood profile of CU *****")
# Profiling over CU
xU,yU,rU = m.mnprofile('CU', size=250, bound=(0., 3), subtract_min=True)

print("***** getting the 1-dimensional likelihood profile of CD *****")
# Profiling over CD
xD,yD,rD = m.mnprofile('CD', size=250, bound=(0., 3), subtract_min=True)

print("***** getting the 1-dimensional likelihood profile of CV *****")
# Profiling over CV
xV,yV,rV = m.mnprofile('CV', size=250, bound=(0., 3), subtract_min=True)


######################################################################
#Exporting data to a .txt file
######################################################################

with open('results/CU_1dprofile_36fb-1.txt', 'w') as f:
	for i in range(0,len(xU)):
		f.write('%.5f    ' % xU[i] +'%.5f     '  % yU[i] + '\n')    
f.close()

with open('results/CD_1dprofile_36fb-1.txt', 'w') as f:
	for i in range(0,len(xD)):
		f.write('%.5f    ' % xD[i] +'%.5f     '  % yD[i] + '\n')
f.close()

with open('results/CV_1dprofile_36fb-1.txt', 'w') as f:
	for i in range(0,len(xV)):
		f.write('%.5f    ' % xV[i] +'%.5f     '  % yV[i] + '\n')
f.close()

