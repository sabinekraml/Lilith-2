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


#outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_20_40_500_150_150_I_Hpm_CMS140fb_2HDMc_CDF.out"
#outputplot = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_20_40_500_150_150_I_Hpm_CMS140fb_2HDMc_CDF.pdf"
#outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_a_tb_60_100_500_100_100_I_Hpm_CMS140fb_2HDMc_PDG.out"
#outputfinal = validation_dir+"mHpm500_CDF_finalresults.out"
#outputplot = validation_dir+"mHpm500_CDF_finalresults_10.pdf"
#outputfinal = validation_dir+"mHpm500_PDG_finalresults.out"
#outputplot = validation_dir+"mHpm500_PDG_finalresults_10.pdf"
#outputfinal = validation_dir + "mHmA_STU_mHpm_minuit_a_tb_30_30_1000_200_100_I_Hpm_CMS140fb_2HDMc.out"
#mH_precision = 30
#mA_precision = 30
outputfinal = validation_dir + "mHmA_STU_mHpm_minuit_a_tb_80_80_1000_200_100_I_Hpm_CMS140fb_2HDMc.out"
outputplot = validation_dir + "mHpm1000_CDF_finalresults_10.pdf"
mH_precision = 80
mA_precision = 80


# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
matplotlib.rcParams['axes.linewidth'] = 2

#fig = plt.figure(figsize=(3.5,7))
fig = plt.figure()
ax = fig.add_subplot(111)

plt.tick_params(direction='out', labelsize=14, length=10, width=2)
#plt.xticks(np.linspace(490, 590, 3))

# Getting the data
data = np.genfromtxt(outputfinal)
x = data[:,0]
y = data[:,1]
z = data[:,2]

#z = np.where(z > np.nanmin(z)+10, np.nan, z)

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
sc = ax.scatter(x, y, c=z2, vmax=10, cmap="jet_r")
#cbar = fig.colorbar(sc,fraction=0.25, pad=0.05)
cbar = fig.colorbar(sc)
cbar.set_label("$\chi^2$", fontsize=16, labelpad=-4)
cbar.ax.tick_params(direction='out', labelsize=14, length=10, width=2)
# 4.695  9.488

# Title, labels, color bar...
#ax.set_xlim([500, 650])
#ax.set_ylim([500, 900])
#ax.set_xlim([490, 590])
#ax.set_ylim([400, 900])
#ax.set_xlim([1000, 1100])
#ax.set_ylim([1100, 1200])
ax.set_xlim([800, 1200])
ax.set_ylim([800, 1200])
ax.set_xlabel(r'$m_H$[GeV]',fontsize=16)
ax.set_ylabel(r'$m_A$[GeV]',fontsize=16)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

plt.show()

print("results are stored in", validation_dir, flush=True)
