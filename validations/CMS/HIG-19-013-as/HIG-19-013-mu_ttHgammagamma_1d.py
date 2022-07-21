###############################################################
# Lilith routine for muttHgammagamma validation plots
###############################################################

import sys, os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np 

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
print(lilith_dir)
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

validation_dir = lilith_dir+"validations/CMS/HIG-19-013-as/"
print(validation_dir)
######################################################################
# Parameters
######################################################################

# Exprimental results
exp_input = validation_dir +"thisRun2.list"


# Lilith precision mode
myprecision = "BEST-QCD"

# Output 
output = validation_dir + "HIG-19-013-muttHgammagamma_1d.out"
outputplot = validation_dir +"HIG-19-013-muttHgammagamma_1d.pdf"

# Scan ranges 
muttHgammagamma_min = 0.
muttHgammagamma_max = 5.

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

# Higgs mass to test
mh = 125.



######################################################################
# * usrXMLinput: generate XML user input
# * getL:        -2LogL for a given muttHgammagamma point
######################################################################


def usrXMLinput(mh=125., muttHgammagamma=1., precision="BEST-QCD"):
    """generate XML input from signal strength muttHgammagamma"""
 
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<signalstrengths>	
  <mass>%(mh)s</mass>
  <mu prod="ttH" decay="gammagamma">%(muttHgammagamma)s</mu>
  <mu prod="ttH" decay="gammagamma">%(muttHgammagamma)s</mu>
  <precision>%(precision)s</precision>
</signalstrengths>

</lilithinput>
"""
    myInput = {'mh':mh, 'muttHgammagamma':muttHgammagamma, 'precision':precision}
        
    return myInputTemplate%myInput
######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

# Prepare output
fresults= open(output, 'w') #Ouvre le fichier stocké dans la variable output pour écrire dessus

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

for muttHgammagamma in np.linspace(muttHgammagamma_min, muttHgammagamma_max, grid_subdivisions):
    fresults.write('\n')
    myXML_user_input = usrXMLinput(mh=mh, muttHgammagamma=muttHgammagamma, precision=myprecision)
    lilithcalc.computelikelihood(userinput=myXML_user_input)
    m2logL = lilithcalc.l
    if m2logL < m2logLmin:
    	m2logLmin = m2logL
    	muttHgammagamma_min = muttHgammagamma
    
    fresults.write('%.5f    ' % muttHgammagamma +'%.5f     '  % m2logL + '\n') #le .5f c'est pour afficher 5 décimales
    
fresults.close()

print("***** scan finalized *****")
print("minimum at muttHgammagamma_min, -2logL_min = ", muttHgammagamma_min, m2logLmin)

######################################################################
# Plot routine
######################################################################

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 8 # distance to major tick label along x axis in points
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(8,8) )
ax = fig.add_subplot(111) #1x1 grid, 1st subplot
ax.set_ylim([0, 50])

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]


#print("data:",data)

# Substracting the -2LogL minimum to form Delta(-2LogL)
y2=[]
for y_el in y:
  y2.append(y_el-y.min())





# read data for official plot 
expdata = np.genfromtxt(validation_dir+'HIG-19-013-Llh-1d-Grid_as.txt')
xexp = expdata[:,0]
yexp = expdata[:,1]

#read data for another mu validation 
datab = np.genfromtxt("data_vn1officialmu.txt")
x_b = datab[:,0]
y2_b = datab[:,1]

#Title, labels, legend... 
plt.xlabel(r'$\mu_{(ttH,\gamma\gamma)}$', fontsize=20)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=20)

#plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's V. Gaussian Appx.")
#plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's Poisson Appx")
#plt.plot(x,y2,'-',markersize=2, color = 'red',label="Lilith V.Gaussian Appx. for $\mu=1.38^{+0.36}_{-0.29}$")
plt.plot(x_b,y2_b,'-',markersize=2, color = 'red',label="Lilith V.Gaussian (1) Appx. for official $\mu=1.38^{+0.36}_{-0.29}$")
plt.plot(x,y2,'-',markersize=2, color = 'g',label="Lilith V.Gaussian (1) Appx. for fitted $\mu=1.43^{+0.34}_{-0.29}$")
plt.plot(xexp,yexp,'.',markersize=3, color = 'blue',label="CMS official values")
plt.legend(loc='upper right', fontsize=12)

#plt.title("$\mu$ from CMS-HIG-19-013 (VGaussian)")
#plt.title("$\mu$ from CMS-HIG-19-013 (Poisson)")
plt.title("Gaussian Variable Appx.(1) for official and fitted $\mu$ from CMS-HIG-19-013")
fig.set_tight_layout(True)

fig.savefig('HIG-19-013-muttHgammagamma_1d_VN1.pdf')
#fig.savefig('mutHgammagamma-Poisson_1D.pdf')
#fig.savefig('mutHgammagamma-grid_1D.pdf')
plt.show()
