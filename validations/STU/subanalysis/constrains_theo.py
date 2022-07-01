import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import subprocess
import time
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor as executor 
from multiprocessing import Pool

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
cba = 0
tb = 1.5

# Precisions
mH_precision = 100
mA_precision = 100

# Lists
mHlist = []
for mH in np.linspace(mH_min, mH_max, mH_precision):
	mHlist.append(mH)

iterationlist = []
for i in range(mH_precision):
	iterationlist.append(i)

# Output files
output = []

# Output files
for i in range(mH_precision):
	output.append(validation_dir+"multiprocessing/constrains_" + str(i) + ".out")

if yukawatype == 1:
	outputfinal = validation_dir+"constrains_theo" + "_" + "mHpmismH" + "_I" + ".out"
	outputplot = validation_dir+"constrains_theo" + "_" + "mHpmismH" + "_I" + ".pdf" 

if yukawatype == 2:
	outputfinal = validation_dir+"constrains_theo" + "_" + "mHpmismH" + "_II" + ".out"
	outputplot = validation_dir+"constrains_theo" + "_" + "mHpmismH" + "_II" + ".pdf" 

######################################################################
# Definition
######################################################################

def func(iteration):
	# Prepare output
	fresults = open(output[iteration], 'w')
	i=0
	mH = mHlist[iteration]

	for mA in np.linspace(mA_min, mA_max, mA_precision):
		if i%(mA_precision/10)==0 and iteration == 0:
			print("mA = ", mA, flush=True)
		i+=1
		mHpm = mH
		b = np.arctan(tb)
		m122 = np.sin(b)**2 * (mH**2/tb)

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(1), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

		if m122>999999:
			if p1.stdout[765] != " ":
				Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[972]), int(p1.stdout[997]), int(p1.stdout[1022])
			else:
				Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
		else:
			Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])


		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(-1), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

		if m122>999999:
			if p1.stdout[765] != " ":
				Treelevelunitarity1, Perturbativity1, Stability1 = int(p1.stdout[972]), int(p1.stdout[997]), int(p1.stdout[1022])
			else:
				Treelevelunitarity1, Perturbativity1, Stability1 = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
		else:
			Treelevelunitarity1, Perturbativity1, Stability1 = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
					
		if ( Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1 ) or ( Treelevelunitarity1 == 1 and Perturbativity1 == 1 and Stability1 == 1 ) :
			fresults.write('%.2f    '%mH + '%.2f    '%mA + '1    ' + '\n')
		else:
			fresults.write('%.2f    '%mH + '%.2f    '%mA + 'nan    ' + '\n')

######################################################################
# Multiprocessing
######################################################################

#if __name__ == '__main__':
#	pool = Pool()
#	pool.map(func, iterationlist)

#fresultsfinal = open(outputfinal, 'w')
#for i in iterationlist:
#	fresults = open(output[i])
#	content = fresults.read()
#	fresultsfinal.write(content+"\n")
#	fresults.close()
#fresultsfinal.close()

######################################################################
# Plot routine
######################################################################

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 6
matplotlib.rcParams['axes.linewidth'] = 2

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.xticks(np.linspace(200, 2000, 7))
plt.yticks(np.linspace(200, 2000, 7))

# Getting the data
data = np.genfromtxt(outputfinal)

x = data[:,0]
y = data[:,1]
z = data[:,2]

# Interpolation
xi = np.linspace(x.min(), x.max(), mH_precision)
yi = np.linspace(y.min(), y.max(), mA_precision)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z, (X, Y), method="linear")

# Plotting
#sc = ax.contourf(xi,yi,Z,[10**(-10),2],colors=['blue'])
sc = ax.scatter(x, y, c=z, s=30)

ax.set_aspect(1)

# Title, labels, color bar...
ax.set_xlabel(r'$m_H = m_{H^{\pm}}$[GeV]',fontsize=16)
ax.set_ylabel(r'$m_A$[GeV]',fontsize=16)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

#plt.show()

print("results are stored in", validation_dir, flush=True)
