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
type = 1

# Scan ranges
mA0 = 1200
mH0 = 1200
mHpm_min = 800
mHpm_max = 1600
m12_min = 0
m12_max = 2000

if type == 1:
  cba_min = -0.25
  cba_max = 0.25
  tb_min = 0.1
  tb_max = 10
if type == 2:
  cba_min = -0.05
  cba_max = 0.05
  tb_min = 0.1
  tb_max = 10

# Precisions
cba_precision = 40
tb_precision = 40
mHpm_precision = 160
m12_precision = 200

cbalist = []
for cba in np.linspace(cba_min, cba_max, cba_precision):
	cbalist.append(cba)

iterationlist = []
for i in range(cba_precision):
	iterationlist.append(i)

# Output files
output = []
outputplot = []
if type == 1:
	for i in range(cba_precision):
		output.append(validation_dir+"multiprocessing/constrains" + "_" + str(i) + ".out")
		
	outputfinal = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(mHpm_precision) + "_" + str(m12_precision) + "_" + "I" + ".out"
	outputplot = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(mHpm_precision) + "_" + str(m12_precision) + "_" + "I" + ".pdf"

if type == 2:
	for i in range(cba_precision):
		output.append(validation_dir+"multiprocessing/constrains" + "_" + str(i) + ".out")
		
	outputfinal = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(mHpm_precision) + "_" + str(m12_precision) + "_" + "II" + ".out"
	outputplot = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(mHpm_precision) + "_" + str(m12_precision) + "_" + "II" + ".pdf"


if type == 2:

	for i in range(cba_precision):
		output.append(validation_dir+"multiprocessing/constrains" + "_" + str(i) + ".out")
		
	outputplot = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(m12_precision) + "_" + "II" + ".pdf"
	outputfinal = validation_dir+"m12/constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm0) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + str(m12_precision) + "_" + "II" + ".out"


######################################################################
# Definition
######################################################################

def func(iteration):

	# Prepare output
	fresults = open(output[iteration], 'w')
	i=0
	cba = cbalist[iteration]

	for tb in np.linspace(tb_min, tb_max, tb_precision):
		if i==1 and iteration==0:
			print("tb = ", tb, flush=True)		
		if i%(tb_precision/10)==0 and iteration==0:
			print("tb = ", tb, flush=True)
		i+=1
		for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
			for m12 in np.linspace(m12_min, m12_max, m12_precision):
				cons = False
				m12_cons = 0
				sba = np.sqrt(1-cba**2)
				p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH0), str(mA0), str(mHpm), str(sba), '0.00000', '0.00000', str(m12), str(tb), str(type)], capture_output=True, text=True)
				Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])

				if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:
					cons = True
					m12_cons = m12
			
			if cons:					
				fresults.write('%.3f    '%cba + '%.2f    '%tb + '%.2f    '%mHpm  + '1    ' + '%.2f    '%m12_cons + '\n')
			else:
				fresults.write('%.3f    '%cba + '%.2f    '%tb + '%.2f    '%mHpm  + 'nan    ' + '%.2f    '%m12_cons + '\n')

	if iteration==0:
		print("tb = ", tb_max, flush=True)
	fresults.close()

######################################################################
# Multiprocessing
######################################################################

start = time.perf_counter()

if __name__ == '__main__':
	pool = Pool()
	pool.map(func, iterationlist)

stop = time.perf_counter()

print("***** scan finalized *****", flush=True)

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
ax.set_xlabel(r'$\cos(\beta - \alpha)$[GeV]',fontsize=10)
ax.set_ylabel(r'$\tan(\beta)$[GeV]',fontsize=10)
ax.set_zlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

#plt.show()

print("results are stored in", validation_dir, flush=True)

#4142216
