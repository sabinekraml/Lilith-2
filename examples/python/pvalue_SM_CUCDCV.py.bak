######################################################################
#
# Lilith routine example
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
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
#myexpinput = "data/latest.list"
myexpinput = "data/latestRun2.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

print "\nTask: Compute the SM and (CU,CD,CV) model p-values"
print "Experimental input:", myexpinput

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

print "\n***** evaluating SM -2LogL *****"
# Evaluate Likelihood at the SM point
SM_minus2logL = getL_CUCDCV(1., 1., 1.)

print "-2LogL(SM) =", SM_minus2logL

print "\n***** performing (CU,CD,CV) model fit *****"
# Initialize the fit; parameter starting values and limits
m = Minuit(getL_CUCDCV, CU=0.9, limit_CU=(0,3), CD=0.9, limit_CD=(0,3), CV=0.9, limit_CV=(0,3), print_level=0, errordef=1, error_CU=1, error_CD=1, error_CV=1)

# Fit the model
m.migrad()
bestfit_CUCDCV_minus2logL = m.fval

# Display parameter values at the best-fit point
print "Best-fit point: "
print " CU =", m.values["CU"], "\n CD =", m.values["CD"], "\n CV =", m.values["CV"] 
print "-2LogL =", bestfit_CUCDCV_minus2logL


print "\n***** statistics *****"

ndf = lilithcalc.exp_ndf
print "Number of degrees of freedom:", ndf

pval_SM = 1 - stats.chi2.cdf(SM_minus2logL, ndf)
pval_CUCDCV = 1 - stats.chi2.cdf(bestfit_CUCDCV_minus2logL, ndf-3)

print "p-value SM:", pval_SM
print "p-value CUCDCV:", pval_CUCDCV, "\n"

