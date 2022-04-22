##################################################################
#
# Lilith routine for (muF, muV) validation plots
#
##################################################################

import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/CMS/HIG-19-006/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

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
#if (not os.path.exists("results")):
#    os.mkdir("results")
output = validation_dir+"HIG-19-006-mumu_2d.out"
outputplot = validation_dir+"HIG-19-006-mumu_2d.pdf"
#output = validation_dir+"HIG-19-006-mumu_2d_cor024.out"
#outputplot = validation_dir+"HIG-19-006-mumu_2d_cor024.pdf"

# Scan ranges 
muf_min = -1.5
muf_max = 2.5
muv_min = 0
muv_max = 5

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100

######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125.38, muf=1, muv=1, precision="BEST-QCD"):
    """generate XML input from signal strengths muf, muv"""
    
    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<signalstrengths>
  <mass>%(mass)s</mass>

  <mu prod="VVH" decay="mumu">%(muv)s</mu>
  <mu prod="ggH" decay="mumu">%(muf)s</mu>
  <mu prod="ttH" decay="mumu">%(muf)s</mu>

  <precision>%(precision)s</precision>
</signalstrengths>

</lilithinput>
"""
    myInput = {'mass':mass, 'muf':muf, 'muv':muv, 'precision':precision}
        
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

for muf in np.linspace(muf_min, muf_max, grid_subdivisions):
    fresults.write('\n')
    for muv in np.linspace(muv_min, muv_max, grid_subdivisions):
        myXML_user_input = usrXMLinput(hmass, muf=muf, muv=muv, precision=my_precision)
        lilithcalc.computelikelihood(userinput=myXML_user_input)
        m2logL = lilithcalc.l
        if m2logL < m2logLmin:
            m2logLmin = m2logL
            mufmin = muf
            muvmin = muv
        fresults.write('%.5f    ' % muf + '%.5f    ' % muv + '%.5f     ' % m2logL + '\n')

fresults.close()

print("***** scan finalized *****")
print("minimum at muf, muv, -2logL_min = ", mufmin, muvmin, m2logLmin)


######################################################################
# Plot routine
######################################################################


print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure()
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

ax.set_aspect((muf_max-muf_min)/(muv_max-muv_min))


# read data for official 68% and 95% CL contours & plot + best data fit point
expdata = np.genfromtxt('validations/CMS/HIG-19-006/HIG-19-006-mumu-Grid.txt')
xExp68 = expdata[1:110,0]
yExp68 = expdata[1:110,1]
plt.plot(xExp68,yExp68,'--',markersize=3, color = '#ff0800', label="CMS official 68% CL")
xExp95 = expdata[110:,0]
yExp95 = expdata[110:,1]
plt.plot(xExp95,yExp95,'--',markersize=3, color = '#ff7b00', label="CMS official 95% CL")
xExpbf = expdata[0,0]
yExpbf = expdata[0,1]
plt.plot(xExpbf,yExpbf,'D',markersize=4, color = '#6cc7e3', label="CMS official best fit")


# best Lilith fit point
plt.plot([mufmin],[muvmin], '*', markersize=8, color = 'black', label = 'Lilith best fit')


# Standard Model 
plt.plot([1], [1], '+',markersize=8, color = '#eed8d7', label="SM prediction")
plt.legend(loc='upper left')


# Title, labels, color bar...
plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
plt.xlabel(r'$\mu$(ggH,ttH)',fontsize=18)
plt.ylabel(r'$\mu$(VBF,VH)',fontsize=18)
#plt.text(0.1, 0.15, r'Data from CMS-HIG-19-006', fontsize=9.5)
#plt.text(0.1, 0.05, r'$\mu_{ggH}$ and $\mu_{top}$ from Fig. 11b', fontsize=9.5)
#plt.text(0.1, -0.05, r'$\mu_{VBF}$ and $\mu_{VH}$ combined in $\mu_{VVH}$ from Fig. 12', fontsize=9.5)
#plt.text(0.1, -0.15, r'$\rho$(ggH,VVH) adjusted to -0.2 from Fig. 15a', fontsize=9.5)


fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("***** done *****")
print("results are stored in", validation_dir)

