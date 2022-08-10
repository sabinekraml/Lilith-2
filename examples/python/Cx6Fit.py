##################################################################
#
# Lilith routine example for (Cx, Cc)  plots
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
import time


lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = "data/latestRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/Fit6parameters.out"
#outputplot = "results/CxCc_2d.pdf"

# Scan ranges
Cc = 1.0
CZ_min = 1.1
CZ_max = 1.2
CW_min = 1.0
CW_max = 1.1
Cb_min = 0.9
Cb_max = 1.1
Ct_min = 0.9
Ct_max = 1.1
Ctau_min = 1
Ctau_max = 1.1
Cmu_min = 0
Cmu_max = 1.5

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 10

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CZ=1, CW=1, Cb=1, Ct=1, Ctau=1, Cc=1, Cmu=1, precision="BEST-QCD"):
    """generate XML input from reduced couplings Cs"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="tt">%(Ct)s</C>
  <C to="bb">%(Cb)s</C>
  <C to="cc">%(Cc)s</C>
  <C to="tautau">%(Ctau)s</C>
  <C to="mumu">%(Cmu)s</C>
  <C to="ZZ">%(CZ)s</C>
  <C to="WW">%(CW)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mass':mass, 'CZ':CZ, 'CW':CW, 'Cb':Cb, 'Ct':Ct, 'Ctau':Ctau, 'Cc':Cc, 'Cmu':Cmu, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

start = time.time()

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

for CZ in np.linspace(CZ_min, CZ_max, grid_subdivisions):
    fresults.write('\n')
    for CW in np.linspace(CW_min, CW_max, grid_subdivisions):
    	fresults.write('\n')
    	for Cb in np.linspace(Cb_min, Cb_max, grid_subdivisions):
    	    fresults.write('\n')
    	    for Ct in np.linspace(Ct_min, Ct_max, grid_subdivisions):
#    	        fresults.write('\n')
    	        for Ctau in np.linspace(Ctau_min, Ctau_max, grid_subdivisions):
#    	            fresults.write('\n')
    	            for Cmu in np.linspace(Cmu_min, Cmu_max, grid_subdivisions):
    	                myXML_user_input = usrXMLinput(hmass, CZ=CZ, CW=CW, Cb=Cb, Ct=Ct, Ctau=Ctau, Cc=Cc, Cmu=Cmu, precision=my_precision)
    	                lilithcalc.computelikelihood(userinput=myXML_user_input)
    	                m2logL = lilithcalc.l
    	                if m2logL < m2logLmin:
    	                    m2logLmin = m2logL
    	                    CZmin = CZ
    	                    CWmin = CW
    	                    Cbmin = Cb
    	                    Ctmin = Ct
    	                    Ctaumin = Ctau
    	                    Cmumin = Cmu
    	    fresults.write('%.5f    '%CZ +'%.5f    '%CW + '%.5f    '%Cb + '%.5f    '%Ct + '%.5f    '%Ctau + '%.5f    '%Cmu + '%.5f     '%m2logL + '\n')
fresults.close()

print("***** scan finalized *****")
print("minimum -2logL_min = ", m2logLmin)
print("at CZ = ",CZmin)
print("   CW = ",CWmin)
print("   Cb = ",Cbmin)
print("   Ct = ",Ctmin)
print(" Ctau = ",Ctaumin)
print("  Cmu = ",Cmumin)
#print("minimum at CZ, CW, Cb, Ct, Ctau, Cc, -2logL_min = ", CZmin, CWmin, Cbmin, Ctmin, Ctaumin, Ccmin, m2logLmin)

end = time.time()
print(" Calculation time = ",end - start)

