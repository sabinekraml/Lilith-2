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

#outputfinal = validation_dir+"constrains_50_50_50_20_20_I.out"
#outputplot = validation_dir+"constrains_50_50_50_20_20_I.pdf"
outputfinal = validation_dir+"constrains_250_300_10_50_50_I.out"
outputplot = validation_dir+"constrains_250_300_10_50_50_I.pdf"

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
ax.set_xlim(600, 700)
ax.set_ylim(-0.25, 0.25)
ax.set_zlim(0, 10)
#cbar = fig.colorbar(sc)

# Title, labels, color bar...
#ax.set_xlabel(r'$m_H$[GeV]',fontsize=10)
#ax.set_ylabel(r'$m_A$[GeV]',fontsize=10)
#ax.set_zlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=10)
ax.set_xlabel(r'$m_{H^{\pm}}$[GeV]',fontsize=10)
ax.set_ylabel(r'$\cos(\beta - \alpha)$[GeV]',fontsize=10)
ax.set_zlabel(r'$\tan(\beta)$[GeV]',fontsize=10)

# Saving figure (.pdf)
fig.savefig(outputplot)

plt.show()

print("results are stored in", validation_dir, flush=True)
