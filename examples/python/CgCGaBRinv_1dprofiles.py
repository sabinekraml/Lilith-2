###############################################################
#
# Lilith routine example
# To put in Lilith-2.X/examples/python/ folder 
# To execute from /Lilith-2.X root folder
#
# Constraints on reduced couplings CGa and Cg 
# and on BR into invisible final states;
# all other couplings as in the SM. 
#
# 1-dimensional likelihood profiles are obtained from a
# profile-likelihood analysis
#
# Uses iminuit (fitting) and matplotlib (plotting)
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
#myexpinput = "data/latest.list"
myexpinput = "data/latestRun2.list"

# Lilith precision mode
myprecision = "BEST-QCD"

# Output
if (not os.path.exists("results")):
    os.mkdir("results")
outputplot = "results/CGaCgBRinv_1dprofiles.pdf"

# Higgs mass to test
myhmass = 125.09

verbose=False
timer=False

######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given (CGa, Cg, BRinv) point
######################################################################

def usrXMLinput(mh=125., CGa=1., Cg=1., BRinv=0., precision="BEST-QCD"):
    """generate XML input from reduced couplings
       the gamma-gamma and/or gluon-gluon couplings are 
       evaluated as functions of the fermionic and bosonic
       reduced couplings """
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mh)s</mass>
  
  <C to="ff">1.</C>
  <C to="VV">1.</C>
  <C to="gammagamma">%(CGa)s</C>
  <C to="gg">%(Cg)s</C>
  
  <extraBR>
    <BR to="invisible">%(BRinv)s</BR>
    <BR to="undetected">0.</BR>
  </extraBR>
  
  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""
    myInput = {'mh':mh, 'CGa':CGa, 'Cg':Cg, 'BRinv':BRinv, 'precision':precision}
        
    return myInputTemplate%myInput



def getL(CGa, Cg, BRinv):
    myXML_user_input = usrXMLinput(mh=myhmass, CGa=CGa, Cg=Cg, BRinv=BRinv, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    return lilithcalc.l


######################################################################
# Calculations
######################################################################

print "***** initializating *****"

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose, timer)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

# Initialize the fit; parameter starting values and limits
m = Minuit(getL, CGa=1, limit_CGa=(0,3), Cg=1, limit_Cg=(0,3), BRinv=0.2, limit_BRinv=(0,0.9), errordef=1, error_CGa=0.1, error_Cg=0.1, error_BRinv=0.1)

print "\n***** performing model fit with iminuit *****"

# Minimization and error estimation
m.migrad()
m.minos()

print "\n***** fit summary *****"

print "\nbest-fit point:", m.values 
print "\nHesse errors:", m.errors
print "\nMinos errors:"
for key, value in m.merrors.items():
    print key, value

print "\nCorrelation matrix:\n", m.matrix(correlation=True)
print "\nCovariance matrix:\n", m.matrix()

# Display parameter values at the best-fit point
#print "Best-fit point: -2LogL =", m.fval
#print "Cgamma =", m.values["CGa"], ", Cgluon =", m.values["Cg"], ", BRinv =", m.values["BRinv"],"\n"

print "\n***** getting the 1d likelihood profiles *****"
# Profiling for CGa
xGa,yGa,rGa = m.mnprofile('CGa', bins=300, bound=(0, 2), subtract_min=True)

# Profiling for Ca
xg,yg,rg = m.mnprofile('Cg', bins=300, bound=(0, 2), subtract_min=True)

# Profiling for BRinv
xBR,yBR,rBR = m.mnprofile('BRinv', bins=300, bound=(0., 0.5), subtract_min=True)


######################################################################
# Plot
######################################################################

print "***** plotting 1d profiles *****"

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)

plt.minorticks_on()
plt.tick_params(labelsize=24, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xGa,yGa,color="b",ls='-',linewidth=2.5,label=r'$C_\gamma$')
plt.plot(xg,yg,color="r",ls='--',linewidth=2.5,label=r'$C_g$')
plt.ylim([0,8])
plt.xlim([0.8,1.5])
plt.xlabel(r'reduced coupling', fontsize=28)
plt.ylabel(r'$-2\,\log\,L/L_0$', fontsize=28)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.legend(loc='upper right', fontsize=24)

plt.title("   Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

ax = fig.add_subplot(122)

plt.minorticks_on()
plt.tick_params(labelsize=24, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xBR,yBR,color="k",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.,0.3])
plt.xlabel(r'BR(invisible)', fontsize=28)
plt.ylabel(r'$-2\,\log\,L/L_0$', fontsize=28)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title("   Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

fig.set_tight_layout(True)
fig.savefig(outputplot)

print  "result see " + lilith_dir + "/" + outputplot
print "***** done *****\n"



