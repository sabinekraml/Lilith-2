import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import subprocess
import time
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor as executor 
from multiprocessing import Pool

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
validation_dir = lilith_dir+"validations/STU/finalanalysis/"
calc2HDM_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))+"/2HDMc/2HDMC-1.8.0/"

#outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_500_I_stra0_2HDMc.out"
#outputplot = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_500_I_stra0_2HDMc.pdf"
#mH_precision = 80
#mA_precision = 80

#outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_40_40_500_100_100_I_stra0_2HDMc.out"
#outputplot = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_40_40_500_100_100_I_stra0_2HDMc_plot.pdf"
outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_40_40_1000_100_100_I_stra0_2HDMc.out"
outputplot = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_40_40_1000_100_100_I_stra0_2HDMc_plot.pdf"
mH_precision = 40
mA_precision = 40

outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_10_10_500_40_40_I_stra0_2HDMc.out"
outputplot = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_10_10_500_40_40_I_stra0_2HDMc_plot.pdf"
mH_precision = 10
mA_precision = 10

# Preparing plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Getting the data
data = np.genfromtxt(outputfinal)
x = data[:,0]
y = data[:,1]
z = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-np.nanmin(z))

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), mH_precision)
yi = np.linspace(y.min(), y.max(), mA_precision)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")

# Plotting
#sc = ax.scatter(x, y, c=z2, cmap="jet_r")
sc = ax.scatter(x, y, c=z2, vmin=0, vmax=10, cmap="jet_r")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$\Delta (-2\log L)$", fontsize=10)

ax.set_aspect(1)

# Title, labels, color bar...
ax.set_xlabel(r'$m_H$[GeV]',fontsize=10)
ax.set_ylabel(r'$m_A$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

plt.show()

print("results are stored in", validation_dir, flush=True)
