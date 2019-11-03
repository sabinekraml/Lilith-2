##################################################################
#
#
##################################################################

import sys, os
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from iminuit import Minuit

lilith_dir = os.getcwd()
sys.path.append(lilith_dir)
import lilith


######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.09, CZ=1., CW=1., CD=1., CU=1., Ctau=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

    <C to="ZZ">%(CZ)s</C>
    <C to="WW">%(CW)s</C>
    <C to="tt">%(CU)s</C>
    <C to="cc">%(CU)s</C>
    <C to="bb">%(CD)s</C>
    <C to="tautau">%(Ctau)s</C>

  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CZ':CZ, 'CW':CW, 'CD':CD, 'CU':CU, 'Ctau':Ctau, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# profile likelihood 
######################################################################

def getL(CZ, CW, CD, CU, Ctau):
    myXML_user_input = usrXMLinput(mass=hmass, CZ=CZ, CW=CW, CD=CD, CU=CU, Ctau=Ctau, precision=my_precision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l

######################################################################
# Scan initialization
######################################################################

my_precision = "BEST-QCD"
hmass = 125.09

lilithcalc = lilith.Lilith(verbose=False,timer=False)

print "***** CZ, CW, CD, CU, Ctau model *****\n"

exp_input = "data/validation.list"
lilithcalc.readexpinput(exp_input)

print "*****************************************"
print "Exp. input:", exp_input 
print "*****************************************"


m1 = Minuit(getL, CZ=1., limit_CZ=(0,2), CW=1., limit_CW=(0,2), CD=1., limit_CD=(0,2), CU=1., limit_CU=(0,2), Ctau=1., limit_Ctau=(0,2), 
  print_level=0, errordef=1, error_CZ=0.05, error_CW=0.05, error_CD=0.05, error_CU=0.05, error_Ctau=0.05)
m1.migrad()
m1.minos()
print "\nbest-fit: "
print m1.values, "\n-2LogL =", m1.fval
#print "\nHesse errors:", m.errors
print "\nMinos errors:"
for key, value in m1.merrors.items():
    print key, value


print "***** done *****\n"




