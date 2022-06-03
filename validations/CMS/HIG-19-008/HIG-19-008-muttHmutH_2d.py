##################################################################
#
# Lilith routine for (mutH, muttH) validation plots
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CMS/HIG-19-008/"

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Experimental results
exp_input = validation_dir+"thisRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125.09

# Output files
#output = validation_dir+"HIG-19-008-mumu_2d_fig14.out"
#outputplot = validation_dir+"HIG-19-008-mumu_2d_fig14.pdf"
output = validation_dir+"HIG-19-008-mumu_2d_fig15a.out"
outputplot = validation_dir+"HIG-19-008-mumu_2d_fig15a.pdf"

# Scan ranges 
mutH_min = -8
mutH_max = 19
muttH_min = -1
muttH_max = 3

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.38, mutH=1, muttH=1, precision="BEST-QCD"):
    """generate XML input from signal strengths mutH, muttH"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<signalstrengths>
  <mass>%(mass)s</mass>

  <mu prod="tH" decay="tautau">%(mutH)s</mu>
  <mu prod="ttH" decay="tautau">%(muttH)s</mu>

  <precision>%(precision)s</precision>
</signalstrengths>

</lilithinput>
"""
    myInput = {'mass':mass, 'mutH':mutH, 'muttH':muttH, 'precision':precision}
        
    return myInputTemplate%myInput

######################################################################
# Scan initialization
######################################################################

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

for mutH in np.linspace(mutH_min, mutH_max, grid_subdivisions):
    fresults.write('\n')
    for muttH in np.linspace(muttH_min, muttH_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, mutH=mutH, muttH=muttH, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            mutHmin = mutH
            muttHmin = muttH
        fresults.write('%.5f    ' % mutH + '%.5f    ' % muttH + '%.5f     ' % m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at mutH, muttH, -2logL_min = ", mutHmin, muttHmin, m2logLmin)


######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())


# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")


# Plotting the 68% and 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99],colors=['#ff3300','#ffa500'])#, \
              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

ax.set_aspect((mutH_max-mutH_min)/(muttH_max-muttH_min))


# read data for official 68% and 95% CL contours & plot + best data fit point
expdata = np.genfromtxt(validation_dir+'HIG-19-008-muttHmutH-Grid.txt')
xExp68 = expdata[1:67,0]
yExp68 = expdata[1:67,1]
plt.plot(xExp68,yExp68,'--',markersize=3, color = '#ff0800', label="CMS official 68% CL")
xExp95 = expdata[67:,0]
yExp95 = expdata[67:,1]
plt.plot(xExp95,yExp95,'--',markersize=3, color = '#ff7b00', label="CMS official 95% CL")
xExpbf = expdata[0,0]
yExpbf = expdata[0,1]
plt.plot(xExpbf,yExpbf,'D',markersize=4, color = '#6cc7e3', label="CMS official best fit")


# best Lilith fit point
plt.plot([mutHmin],[muttHmin], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\mu$(tH)',fontsize=18)
plt.ylabel(r'$\mu$(ttH)',fontsize=18)
plt.text(-6.5, -0.6, r'Data from CMS-HIG-19-008', fontsize=9.5)
#plt.text(-6.5, -0.8, r'$\mu$ from Fig. 14 and $\rho$ fitted from 95% CL in Fig. 15a', fontsize=9.5)
plt.text(-6.5, -0.8, r'$\mu$ and $\rho$ fitted from 95% CL in Fig. 15a', fontsize=9.5)


fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", validation_dir)

