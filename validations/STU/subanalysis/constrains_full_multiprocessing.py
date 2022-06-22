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
mHpm_min = 200
mHpm_max = 2000

#if yukawatype == 1:
#  cba_min = -0.25
#  cba_max = 0.25
#  tb_min = 0.1
#  tb_max = 10
#if yukawatype == 2:
#  cba_min = -0.05
#  cba_max = 0.05
#  tb_min = 0.1
#  tb_max = 10

a_min = 0
a_max = np.pi/2
tb_min = 0.5
tb_max = 10

# Precisions
mH_precision = 40
mA_precision = 40
mHpm_precision = 40
#cba_precision = 20
a_precision = 20
tb_precision = 20
#mH_precision = 2
#mA_precision = 2
#mHpm_precision = 2
#cba_precision = 2
#tb_precision = 2

# Lists
mHlist = []
for mH in np.linspace(mH_min, mH_max, mH_precision):
	mHlist.append(mH)

iterationlist = []
for i in range(mH_precision):
	iterationlist.append(i)

# Output files
output = []
outputplot = []
if yukawatype == 1:
	for i in range(mH_precision):
		output.append(validation_dir+"multiprocessing/constrains_" + str(i) + ".out")
#		outputplot = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "I" + "_" + str(i) + ".pdf"

#	outputfinal = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "I" + ".out"
		outputplot = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_" + "I" + "_" + str(i) + ".pdf"

	outputfinal = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_" + "I" + ".out"


if yukawatype == 2:
	for i in range(mH_precision):
		output = validation_dir+"multiprocessing/constrains_" + str(i) + ".out"
#		outputplot = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "II" + "_" + str(i) + ".pdf"

#	outputfinal = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "II" + ".out"
		outputplot = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_" + "II" + "_" + str(i) + ".pdf"

	outputfinal = validation_dir+"constrains" + "_" + str(mA_precision) + "_" + str(mH_precision) + "_" + str(mHpm_precision) + "_" + str(a_precision) + "_" + str(tb_precision) + "_" + "II" + ".out"

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
		for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
#			for cba in np.linspace(cba_min, cba_max, cba_precision):
			cons = False
			for a in np.linspace(a_min, a_max, a_precision):
#				print("cba = ", cba)
#				print("cons = ", cons)
				for tb in np.linspace(tb_min, tb_max, tb_precision):
#					sba = np.sqrt(1-cba**2)
#					m122 = ( np.sin(np.arctan(tb))*sba + cba*np.cos(np.arctan(tb)) )**2 * (mH**2/tb)
					sba = np.sin(np.arctan(tb)-a)
					m122 = np.cos(a)**2*mH**2/tb

					p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(sba), '0.00000', '0.00000', str(m122), str(tb), str(yukawatype)], capture_output=True, text=True)

					if m122>999999:
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[971]), int(p1.stdout[996]), int(p1.stdout[1021])
					else:
						Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
			
					if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1				
						fresults.write('%.2f    '%mH + '%.2f    '%mA + '%.2f    '%mHpm + '1    ' + '%.2f    '%cba + '%.2f    '%tb + '\n')
						cons = True
						break

				if cons:
					break

			if cons=False:
			fresults.write('%.2f    '%mH + '%.2f    '%mA + '%.2f    '%mHpm + 'nan    ' + 'nan    ' + 'nan    ' + '\n')

		if iteration == 0 and mA is not mA_max:
			print("time = ", time.perf_counter()-start, flush=True)
			
	if iteration == 0:
		print("mA = ", mA, flush=True)
		print("time = ", time.perf_counter()-start, flush=True)
		
	print("mH = ", mH, " : done")

	fresults.close()

######################################################################
# Multiprocessing
######################################################################

start = time.perf_counter()

if __name__ == '__main__':
	pool = Pool()
	pool.map(func, iterationlist)

stop = time.perf_counter()

print("time = ", stop-start, flush=True)

fresultsfinal = open(outputfinal, 'w')
for i in iterationlist:
	fresults = open(output[i])
	content = fresults.read()
	fresultsfinal.write(content+"\n")
	fresults.close()
fresultsfinal.close()

######################################################################
# Plot routine
######################################################################

# Preparing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Getting the data
data = np.genfromtxt(outputfinal)
x = data[:,0]
y = data[:,1]
z = data[:,2]
consvalue = data[:,3]

# Plotting
sc = ax.scatter(x, y, z, c=consvalue, s=30)
#cbar = fig.colorbar(sc)

# Title, labels, color bar...
ax.set_xlabel(r'$m_H$[GeV]',fontsize=10)
ax.set_ylabel(r'$m_A$[GeV]',fontsize=10)
ax.set_zlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

#plt.show()

print("results are stored in", validation_dir, flush=True)
