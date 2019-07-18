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
myexpinput = "data/validation.list"
# Lilith precision mode
myprecision = "BEST-QCD"
# Output
plottitle = "validation/CMS/CGaCgBRinv_profiles.pdf"
# Higgs mass to test
myhmass = 125.

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

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose, timer)
# Read experimental data
lilithcalc.readexpinput(myexpinput)

# Initialize the fit; parameter starting values and limits
#m = Minuit(getL, CGa=1, limit_CGa=(0,5), Cg=1, limit_Cg=(0,5), BRinv=0.2, limit_BRinv=(0,1), print_level=0, errordef=0.5, error_CGa=0.5, error_Cg=0.5, error_BRinv=0.1)
m = Minuit(getL, CGa=1, limit_CGa=(0,2), Cg=1, limit_Cg=(0,2), BRinv=0.2, limit_BRinv=(0,0.6), errordef=1, error_CGa=0.1, error_Cg=0.1, error_BRinv=0.1)

# Minimization and error estimation
m.migrad()
m.minos()

# Display parameter values at the best-fit point
print "Best-fit point: -2LogL =", m.fval
print "Cgamma =", m.values["CGa"], ", Cgluon =", m.values["Cg"], ", extraBR =", m.values["BRinv"],"\n"

#print(m.get_param_states(),"\n\n")
#print(m.hesse(),"\n\n")
print "Correlation matrix:\n", m.matrix(correlation=True),"\n"

print "***** getting the 1-dimensional likelihood profiles *****"
# Profiling over CGa
xGa,yGa,rGa = m.mnprofile('CGa', bins=300, bound=(0.1, 2), subtract_min=True)

# Profiling over Ca
xg,yg,rg = m.mnprofile('Cg', bins=300, bound=(0.1, 2), subtract_min=True)

# Profiling over extraBR
xBR,yBR,rBR = m.mnprofile('BRinv', bins=300, bound=(0., 0.6), subtract_min=True)


######################################################################
# Plot
######################################################################

print "***** plotting 1d profiles *****\n"

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)

plt.minorticks_on()
plt.tick_params(labelsize=24, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xGa,yGa,color="b",linewidth=2.5,label=r'$C_\gamma$')
plt.plot(xg,yg,color="r",linewidth=2.5,label=r'$C_g$')
plt.ylim([0,8])
plt.xlim([0.8,1.5])
plt.xlabel(r'reduced coupling', fontsize=28)
plt.ylabel(r'$-2\,\log\,L/L_0$', fontsize=28)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.legend(loc='upper right', fontsize=24)

plt.title(" Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=21, ha="left")

ax = fig.add_subplot(122)

plt.minorticks_on()
plt.tick_params(labelsize=24, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xBR,yBR,color="k",linewidth=2.5)
plt.ylim([0,8])
plt.xlim([0.,0.5])
plt.xlabel(r'BR(invisible)', fontsize=28)
plt.ylabel(r'$-2\,\log\,L/L_0$', fontsize=28)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')

plt.title(" Lilith "+str(lilith.__version__)+", DB "+str(lilithcalc.dbversion), fontsize=20, ha="left")

#plt.tight_layout()
#plt.savefig(plottitle, bbox_inches='tight')
fig.set_tight_layout(True)
fig.savefig("validation/CMS/HIG-17-031-CgCGaBRinv-1dprofiles-v2.pdf")



