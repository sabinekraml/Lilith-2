###############################################################
#
# Lilith routine example
# To execute from /Lilith-1.X root folder
#
# Compute the SM and (CU,CD,CV) model p-values
#
# Use the library iminuit to perform minimization
# Use the library scipy to compute p-values
#
###############################################################

import sys, os
from iminuit import Minuit
from math import floor
from scipy import stats

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith


######################################################################
# Parameters
######################################################################

# Choose which experimental data to fit
myexpinput = "data/latest.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Higgs mass to test
myhmass = 125.

verbose=False
timer=False


######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given (CU,CD,CV)
######################################################################

def usrXMLinput(mass=125., CU=1., CD=1., CV=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings
       CU, CD, CV"""
    
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
    myInput = {'mass':mass, 'CU':CU, 'CD':CD, 'CV':CV, 'precision':precision}
        
    return myInputTemplate%myInput



def getL_CUCDCV(CU, CD, CV):
  """ 
      Returns -2Log(L) for a model in which up-type fermions,
      down-type fermions and gauge bosons have common coupling 
      CU, CD, CV respectively. The loop-induced gamma-gamma and
      gluon-gluon couplings are derived from CU, CD, CV. 
  """
  myXML_user_input = usrXMLinput(myhmass, CU=CU, CD=CD, CV=CV, precision=myprecision)
  lilithcalc.computelikelihood(userinput=myXML_user_input)
  return lilithcalc.l


######################################################################
# Calculations
######################################################################

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose,timer)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

print "***** evaluating SM -2LogL *****"
# Evaluate Likelihood at the SM point
SM_minus2logL = getL_CUCDCV(1., 1., 1.)

print "***** initializing (CU,CD,CV) model fit *****"
# Initialize the fit; parameter starting values and limits
m = Minuit(getL_CUCDCV, CU=0.9, limit_CU=(0,3), CD=0.9, limit_CD=(0,3), CV=0.9, limit_CV=(0,3), print_level=0, errordef=1, error_CU=1, error_CD=1, error_CV=1)

print "***** performing (CU,CD,CV) model fit *****"
# Fit the model
m.migrad()
# Display parameter values at the best-fit point
print "\nBest-fit point of the (CU,CD,CV) model: "
print "CU:", m.values["CU"], ", CD:", m.values["CD"], ", CV:", m.values["CV"], "\n"

bestfit_CUCDCV_minus2logL = m.fval

ndf = lilithcalc.exp_ndf
print "Number of degrees of liberty :", ndf
print "-2LogL(SM) =", SM_minus2logL
print "-2LogL(best-fit CUCDCV model) =", bestfit_CUCDCV_minus2logL

pval_SM = 1 - stats.chi2.cdf(SM_minus2logL, ndf)
pval_CUCDCV = 1 - stats.chi2.cdf(bestfit_CUCDCV_minus2logL, ndf-3)

print "p-value SM:", pval_SM
print "p-value CUCDCV:", pval_CUCDCV

