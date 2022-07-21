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
output = np.genfromtxt("results/CV_1dprofile_36fb-1.txt")
x = output[:,0]
y = output[:,1]
outputplot = "plots/CV-CUCD_1dprofile_36fb-1.pdf"

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

 

plt.plot(x,y,color="b",linewidth=2.5)
plt.ylim([0,10])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_U$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.text(1.16, 2,"Min. at 1.08 $\pm$ 0.06", fontsize=20)

plt.title("Lilith-2.1, 36 fb$^{-1}$", fontsize=23, ha="center")

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("***** done *****\n")
