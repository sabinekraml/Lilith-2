import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import subprocess
import time

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
validation_dir = lilith_dir+"validations/STU/subanalysis/"
calc2HDM_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))+"/2HDMc/2HDMC-1.8.0/"


######################################################################
# Parameters
######################################################################

# 2HDM type = 1, 2
yukawatype = 1

# Scan ranges
mA_min = 200
mA_max = 2000
mH_min = 200
mH_max = 2000
mHpm = 1000
cba = 0
tb = 1.5

# Precisions
mA_precision = 100
mH_precision = 100

# Output files
if yukawatype == 1:
	output = validation_dir+"constrains_theo" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm) + "_" + str(cba) + "_" + str(tb) + "_I" + ".out"
	outputplot = validation_dir+"constrains_theo"  + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm) + "_" + str(cba) + "_" + str(tb) + "_I" + ".pdf"

if yukawatype == 2:
	output = validation_dir+"constrains_theo" + "_"  + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm) + "_" + str(cba) + "_" + str(tb) + "_II" + ".out"
	outputplot = validation_dir+"constrains_theo" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm) + "_" + str(cba) + "_" + str(tb) + "_II" +".pdf"


######################################################################
# Definition
######################################################################

# Prepare output
fresults = open(output, 'w')

i=0

for mH in np.linspace(mH_min, mH_max, mH_precision):
	print("mH = ", mH, flush=True)
	for mA in np.linspace(mA_min, mA_max, mA_precision):
		print("mA = ", mA)
		b = np.arctan(tb)
		sinba = np.sqrt(1-cba**2)
		m122 = ( np.sin(b)*sinba + cba*np.cos(b) )**2 * (mH**2/tb)

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(sinba), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

		if m122>999999:
			if p1.stdout[765] != " ":
				Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[972]), int(p1.stdout[997]), int(p1.stdout[1022])
			else:
				Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
		else:
			Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
					
		if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:
			fresults.write('%.2f    '%mH + '%.2f    '%mA + '1    ' + '\n')
		else:
			fresults.write('%.2f    '%mH + '%.2f    '%mA + 'nan    ' + '\n')

######################################################################
# Plot routine
######################################################################

# Preparing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,2]
consvalue = data[:,4]

# Plotting
sc = ax.scatter(x, y, z, c=consvalue, s=30)

# Title, labels, color bar...
ax.set_xlabel(r'$m_H$[GeV]',fontsize=10)
ax.set_ylabel(r'$m_A$[GeV]',fontsize=10)


# Saving figure (.pdf)
fig.savefig(outputplot)

plt.show()

print("results are stored in", validation_dir, flush=True)
