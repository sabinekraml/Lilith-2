import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import subprocess

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"

validation_dir = lilith_dir+"validations/STU/subanalysis/"

calc2HDM_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))+"/2HDMc/2HDMC-1.8.0/"

# Values
mA0 = 1600
mH0 = 1100
mHpm0 = 350
tb0 = 0.3

# Scan ranges
cba_min = -0.25
cba_max = 0.25

# Precision
precision = 100

# 2HDM type
yukawatype = 2

# Output files
output = validation_dir+"2HDMc_cba_" + str(mA0) + "_" + str(mH0) + "-" + str(mHpm0) + "-" + str(tb0) + "-" + str(yukawatype) + ".out"
outputplot = validation_dir+"2HDMc_cba_" + str(mA0) + "_" + str(mH0) + "-" + str(mHpm0) + "-" + str(tb0) + "-" + str(yukawatype) + ".pdf"

fresults = open(output, 'w')

i=1

for cba in np.linspace(cba_min, cba_max, precision):
		if i==1:
				print("cba = ", cba, flush=True)
		if i%10==0:
				print("cba = ", cba, flush=True)
		i+=1
		cons = False
		m12 = np.cos( np.arctan(tb0) - np.arccos(cba) ) * (mH0/np.sqrt(tb0))
		sinba = np.sqrt(1-cba**2)

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH0), str(mA0), str(mHpm0), str(sinba), '0.00000', '0.00000', str(m12), str(tb0), str(yukawatype)], capture_output=True, text=True)
		S, T, U = float(p1.stdout[1056:1068]), float(p1.stdout[1083:1095]), float(p1.stdout[1110:1122])
		Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])

		cons = Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1
		fresults.write('%.5f    '%cba + '%.1f    '%Treelevelunitarity + '%.1f    '%Perturbativity + '%.1f    '%Stability + '%.1f    '%cons + '%.1f    '%S + '%.1f    '%T + '%.1f    '%U + '\n')
