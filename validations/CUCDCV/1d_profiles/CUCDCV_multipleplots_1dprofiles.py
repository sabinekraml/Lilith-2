import sys, os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np 

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

if (not os.path.exists("plots")):
    os.mkdir("plots")

######################################################################
# Plot
######################################################################

# Output 
outputCU = np.genfromtxt("results/CU_1dprofile_140fb-1.txt")
outputCD = np.genfromtxt("results/CD_1dprofile_140fb-1.txt")
outputCV = np.genfromtxt("results/CV_1dprofile_140fb-1.txt")
outputplot = "plots/CUCDCV_1dprofiles_140fb-1.pdf"

xU = outputCU[:,0]
yU = outputCU[:,1]
xD = outputCD[:,0]
yD = outputCD[:,1]
xV = outputCV[:,0]
yV = outputCV[:,1] 

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(221)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

 

plt.plot(xU,yU,color="b",linewidth=2.5)
plt.ylim([0,10])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_U$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.text(1.13, 2,"Min. at 0.96 $\pm$ 0.04", fontsize=20)

plt.title("Lilith-2.1, 140 fb$^{-1}$", fontsize=23, ha="center")

ax = fig.add_subplot(222)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xD,yD,color="b",linewidth=2.5)
plt.ylim([0,10])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_D$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.text(1.13, 2,"Min. at 0.98 $\pm$ 0.08", fontsize=20)

plt.title("Lilith-2.1, 140 fb$^{-1}$", fontsize=23, ha="center")

ax = fig.add_subplot(223)
plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.plot(xV,yV,color="b",linewidth=2.5)
plt.ylim([0,10])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_V$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.text(1.13, 2,"Min. at 1.04 $\pm$ 0.04", fontsize=20)

plt.title("Lilith-2.1, 140 fb$^{-1}$", fontsize=23, ha="center")
fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("***** done *****\n")
