import sys, os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np 

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))+"/"
sys.path.append(lilith_dir)
import lilith

if (not os.path.exists("plots")):
    os.mkdir("plots")

######################################################################
# Plot
######################################################################

# Output 
output36 = np.genfromtxt(lilith_dir+"validations/CUCDCV/1d_profiles/results/CD_1dprofile_36fb-1.txt")
output140 = np.genfromtxt(lilith_dir+"validations/CUCDCV/1d_profiles/results/CD_1dprofile_140fb-1.txt")
x36 = output36[:,0]
y36 = output36[:,1]
x140 = output140[:,0]
y140 = output140[:,1]
outputplot = "plots/CD-CUCV_1dprofile.pdf"

print("***** plotting *****")

matplotlib.rcParams['xtick.major.pad'] = 18
matplotlib.rcParams['ytick.major.pad'] = 18

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

plt.minorticks_on()
plt.tick_params(labelsize=26, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

 

plt.plot(x36,y36,color="b",linewidth=2.5,label="36 fb$^{-1}$ data")
plt.plot(x140,y140,color="r",linewidth=2.5,label="140 fb$^{-1}$ data")
plt.ylim([0,10])
plt.xlim([0.7,1.4])
plt.xlabel(r'$C_D$', fontsize=34)
plt.ylabel(r'$\Delta (-2\log L)$', fontsize=34)
plt.axhline(y=1.,color='k',ls='dashed')
plt.axhline(y=4.,color='k',ls='dashed')
plt.axhline(y=9.,color='k',ls='dashed')
plt.text(0.72, 2,"Min. at 1.06 $\pm$ 0.11", fontsize=20,color="blue")
plt.text(0.72, 1.2,"Min. at 0.98 $\pm$ 0.08", fontsize=20,color="red")

plt.title("Lilith-2.1, 36 and 140 fb$^{-1}$ data comparison", fontsize=23, ha="center")
plt.legend(loc="upper left",fontsize=20)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot, bbox_inches='tight')

print("***** done *****\n")
