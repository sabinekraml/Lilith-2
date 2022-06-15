import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import subprocess

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
validation_dir = lilith_dir+"validations/STU/subanalysis/"
calc2HDM_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))+"/2HDMc/2HDMC-1.8.0/"


######################################################################
# Parameters
######################################################################

# 2HDM type = 1, 2
type = 1

# Scan ranges
mA0 = 250
mH0 = 300
mHpm_min = 600
mHpm_max = 700
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
mHpm_precision = 10
cba_precision = 50
tb_precision = 50

# Output files
if type == 1:
  output = validation_dir+"constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "I" + ".out"
  outputplot = validation_dir+"constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "I" + ".pdf"

if type == 2:
  output = validation_dir+"constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "II" + ".out"
  outputplot = validation_dir+"constrains" + "_" + str(mA0) + "_" + str(mH0) + "_" + str(mHpm_precision) + "_" + str(cba_precision) + "_" + str(tb_precision) + "_" + "II" + ".pdf"


# Prepare output
fresults = open(output, 'w')


i=0

for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
	if i%(mHpm_precision/10)==0:
			print("mHpm = ", mHpm, flush=True)
	i+=1
#	print("mHpm = ", mHpm, flush=True)
	for cba in np.linspace(cba_min, cba_max, cba_precision):
		for tb in np.linspace(tb_min, tb_max, tb_precision):
			m12 = np.cos( np.arctan(tb) - np.arccos(cba) ) * (mH0/np.sqrt(tb))
			sba = np.sqrt(1-cba**2)
			p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH0), str(mA0), str(mHpm), str(sba), '0.00000', '0.00000', str(m12), str(tb), str(type)], capture_output=True, text=True)
			Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
				
			if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:				
				fresults.write('%.2f    '%mHpm + '%.2f    '%cba + '%.2f    '%tb + '1    ' + '\n')
			else:
					fresults.write('%.2f    '%mHpm + '%.2f    '%cba + '%.2f    '%tb+ 'nan    ' + '\n')

print("mHpm = ", mHpm, flush=True)
fresults.close()

print("***** scan finalized *****", flush=True)

######################################################################
# Plot routine
######################################################################

## Preparing plot
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

## Getting the data
#data = np.genfromtxt(output)
#x = data[:,0]
#y = data[:,1]
#z = data[:,2]
#consvalue = data[:,3]

## Plotting
#sc = ax.scatter(x, y, z, c=consvalue)
##cbar = fig.colorbar(sc)

## Title, labels, color bar...
#ax.set_xlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=10)
#ax.set_ylabel(r'$\cos(\beta - \alpha)$[GeV]',fontsize=10)
#ax.set_zlabel(r'$\tan(\beta)$[GeV]',fontsize=10)

## Saving figure (.pdf)
#fig.savefig(outputplot)

##plt.show()

#print("results are stored in", validation_dir)

print("***** done *****", flush=True)
