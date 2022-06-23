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

angletype = "a"

# Scan ranges
mass_min = 200
mass_max = 2000
cba_min = 0
cba_max = 1
a_min = -np.pi/2
a_max = np.pi/2
tb_min = 0.5
tb_max = 10

# Precisions
mass_precision = 20
cba_precision = 20
a_precision = 20
tb_precision = 20

# Output files
if angletype=="a":
	if yukawatype == 1:
		output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_I_" + "a_range" + ".out"
		outputplot = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_I_" + "a_range" + ".pdf"

	if yukawatype == 2:
		output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_II_" + "a_range" + ".out"
		outputplot = validation_dir+"constrains_minmial" + "_" + str(mass_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_II_" + "a_range" + ".pdf"

if angletype=="cba":
	if yukawatype == 1:
		output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_I_" + "cba" + ".out"
		outputplot = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_I_" + "cba" + ".pdf"

	if yukawatype == 2:
		output = validation_dir+"constrains_minimal" + "_" + str(mass_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_II_" + "cba" + ".out"
		outputplot = validation_dir+"constrains_minmial" + "_" + str(mass_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_II_" + "cba" + ".pdf"

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
	if angletype=="a":
		for a in np.linspace(a_min, a_max, a_precision):
			for tb in np.linspace(tb_min, tb_max, tb_precision):
				b = np.arctan(tb)
				sinba = np.sin(b-a)
				m122 = np.cos(a)**2*mass**2/tb

				p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mass), str(mass), str(mass), str(sinba), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

				if m122>999999:
#					print(p1.stdout)
					if p1.stdout[765] != " ":
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[972]), int(p1.stdout[997]), int(p1.stdout[1022])
					else:
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
				else:
					Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
					
				if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:				
					fresults.write('%.2f    '%mass + '%.4f    '%a + '%.4f    '%tb + '1    ' + '%.2f    '%m122 + '\n')
				else:
					fresults.write('%.2f    '%mass + '%.4f    '%a + '%.4f    '%tb + 'nan    ' + 'nan    ' + '\n')
	else:
		for cba in np.linspace(cba_min, cba_max, cba_precision):
			for tb in np.linspace(tb_min, tb_max, tb_precision):
				b = np.arctan(tb)
				sinba = np.sqrt(1-cba**2)
				m122 = ( np.sin(b)*sinba + cba*np.cos(b) )**2 * (mass**2/tb)

				p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mass), str(mass), str(mass), str(sinba), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

				if m122>999999:
					if p1.stdout[765] != " ":
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[972]), int(p1.stdout[997]), int(p1.stdout[1022])
					else:
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
				else:
					Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
					
				if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:				
					fresults.write('%.2f    '%mass + '%.4f    '%cba + '%.4f    '%tb + '1    ' + '%.2f    '%m122 + '\n')
				else:
					fresults.write('%.2f    '%mass + '%.4f    '%cba + '%.4f    '%tb + 'nan    ' + 'nan    ' + '\n')


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
ax.set_xlabel(r'$m_H = m_A = m_{{H^{{\pm}}}}$[GeV]',fontsize=10)
if angletype=="a":
	ax.set_ylabel(r'$\alpha$[rad]',fontsize=10)
if angletype=="cba":
	ax.set_ylabel(r'$\cos(\beta - \alpha)$[rad]',fontsize=10)
ax.set_zlabel(r'$\tan(\beta)$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

plt.show()

print("results are stored in", validation_dir, flush=True)
