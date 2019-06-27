###############################################################
#
# Lilith routine example
#
# To execute from /Lilith-1.X root folder
#
# Constraints on reduced couplings
# in the (CV, CF) model. CV: common coupling to EW gauge bosons,
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
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith


######################################################################
# Parameters
######################################################################

# Exprimental results
#myexpinput = "data/ATLAS_36fb_full.list"
#myexpinput = "data/ATLAS_best.list"
myexpinput = "data/test.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Output
plottitle = "results/CVCF_1dprofile_ATLAS_1Df_new.pdf"
# Higgs mass to test
myhmass = 125.

verbose=False
timer=False

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

print "***** initializing (CV, CF) model fit *****"
# Initialize the fit; parameter starting values and limits
m = Minuit(getL, CV=1, limit_CV=(0,3), CF=1, limit_CF=(0,3), print_level=0, errordef=1, error_CV=1, error_CF=1)

print "***** performing (CV, CF) model fit *****"
# Fit the model
m.migrad()
# Display parameter values at the best-fit point
print "\nBest-fit point of the (CV, CF) model: "
print "CV =", m.values["CV"], ", CF =", m.values["CF"],"\n"

print "***** getting the 1-dimensional likelihood profile of CV *****"
# Profiling over CF
xV,yV,rV = m.mnprofile('CV', bins=300, bound=(0., 3), subtract_min=True)

print "***** getting the 1-dimensional likelihood profile of CF *****\n"
# Profiling over CV
xF,yF,rF = m.mnprofile('CF', bins=300, bound=(0., 3), subtract_min=True)


######################################################################
# Plot
######################################################################

print "***** plotting *****"

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)

plt.minorticks_on()
plt.tick_params(labelsize=28, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xV,yV,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.8,1.3])
plt.xlabel(r'$C_V$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title(" Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=21, ha="left")

ax = fig.add_subplot(122)

plt.minorticks_on()
plt.tick_params(labelsize=25, length=10, width=1.3)
plt.tick_params(which='minor', length=6, width=0.8)

plt.plot(xF,yF,color="b",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.6,1.501])
plt.xlabel(r'$C_F$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title(" Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=21, ha="left")

plt.tight_layout()
plt.savefig(plottitle, bbox_inches='tight')


