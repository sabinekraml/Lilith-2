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

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)
import iminuit
print("iminuit version:", iminuit.__version__)

######################################################################
# Parameters
######################################################################

# Exprimental results
myexpinput = "../../data/latestRun2.list"

# Lilith precision mode
myprecision = "BEST-QCD"

# Output
outputplot = "CVCF_1dprofiles_minuit2.pdf"

# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

print("\nTask: profile-likelihood analysis in (CF,CV) model")
print("Experimental input:", myexpinput)

######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given (CV, CF) point
######################################################################

def usrXMLinput(mh=125., CV=1., CF=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings
       the gamma-gamma and/or gluon-gluon couplings are 
       evaluated as functions of the fermionic and bosonic
       reduced couplings """
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mh)s</mass>
  
  <C to="tt">%(CF)s</C>
  <C to="bb">%(CF)s</C>
  <C to="tautau">%(CF)s</C>
  <C to="cc">%(CF)s</C>
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
    myInput = {'mh':mh, 'CV':CV, 'CF':CF, 'precision':precision}
        
    return myInputTemplate%myInput



def getL(CV, CF):
    myXML_user_input = usrXMLinput(mh=myhmass, CV=CV, CF=CF, precision=myprecision)
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
CV0 = 1
CF0 = 1

# Initialize the fit; parameter starting values and limits
#least_squares = LeastSquares(CV, CF, )
#m = Minuit(getL, CV=CV0, CF=CF0, print_level=0, errordef=1, error_CV=0.2, error_CF=0.2)
m = Minuit(getL, CV=CV0, CF=CF0)

m.limits = [(0,3), (0,3)]

# Minimization and error estimation
m.migrad()
m.hesse()   # run covariance estimator
m.minos()


print("\nbest-fit point:", m.values) 
print("\nHesse errors:", m.errors)
print("\nMinos errors:")
for key, value in list(m.merrors.items()):
    print(key, value)

#print("\nCorrelation matrix:\n", m.matrix(correlation=True))
#print("\nCovariance matrix:\n", m.matrix())


# Display parameter values at the best-fit point
#print "\nBest-fit point:" 
#print "CV =", m.values["CV"], "\nCF =", m.values["CF"], "\n-2LogL =", m.fval


print("\n***** getting the 1-dimensional likelihood profile of CV *****")
# Profiling over CF
xV,yV,rV = m.mnprofile('CV', size=100, bound=(0., 3), subtract_min=True)

print("***** getting the 1-dimensional likelihood profile of CF *****")
# Profiling over CV
xF,yF,rF = m.mnprofile('CF', size=100, bound=(0., 3), subtract_min=True)


######################################################################
# Plot
######################################################################

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xV,yV,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.8,1.3])
plt.xlabel(r'$C_V$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

ax = fig.add_subplot(122)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xF,yF,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_F$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title("Lilith-"+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

#fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("result see " + validation_dir)
print("***** done *****\n")