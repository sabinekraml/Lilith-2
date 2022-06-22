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
mass_min = 1000
mass_max = 2000
m122_min = 10000
m122_max = 3000000

cba0 = 0
sba0 = 1
tb0 = 1.5

# Precisions
mass_precision = 10
m122_precision = 300

# Output files
if yukawatype == 1:
	output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(m122_precision) + "_" + "I" + ".out"
	outputplot = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(m122_precision) + "_" + "I" + ".pdf"

if yukawatype == 2:
	output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(m122_precision) + "_" + "II" + ".out"
	outputplot = validation_dir+"constrains_minmial" + "_" + str(mass_precision) + "_" + str(m122_precision) + "_" + "II" + ".pdf"

######################################################################
# Definition
######################################################################

# Prepare output
fresults = open(output, 'w')

i=0

for mass in np.linspace(mass_min, mass_max, mass_precision):
	j=0
	if i%(mass_precision/10)==0:
		print("mass = ", mass, flush=True)
	i+=1
	for m122 in np.linspace(m122_min, m122_max, m122_precision):
#		if j%(m122_precision/10)==0:
#			print("m122 = ", m122, flush=True)
#		j+=1
		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mass), str(mass), str(mass), str(-1), '0.00000', '0.00000', str(m122), str(1.5), str(yukawatype)], capture_output=True, text=True)

		if m122>999999:
			Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
		else:
			Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
			
		if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:				
			fresults.write('%.4f    '%mass + '%.2f    '%m122 + '1    ' + '\n')
		else:
			fresults.write('%.4f    '%mass + '%.2f    '%m122 + 'nan    ' + '\n')


######################################################################
# Plot routine
######################################################################

# Preparing plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Getting the data
data = np.genfromtxt(output)
x = data[:,0]
y = data[:,1]
consvalue = data[:,2]

# Plotting
sc = ax.scatter(x, y, c=consvalue, s=30)

# Title, labels, color bar...
ax.set_xlabel(r'$m_H = m_A = m_{{H^{{\pm}}}}$[GeV]',fontsize=10)
ax.set_ylabel(r'$m_{12}$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

#plt.show()

print("results are stored in", validation_dir, flush=True)
