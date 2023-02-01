# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
#lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

# ====================================================================== 
#			declaring input and output  
# ======================================================================
# Experimental results
exp_input = "validations/ATLAS/HIGG-2018-28-SMEFT/latestRun2-stxs.list"

# Sm predictions     
smpred_input = "validations/ATLAS/HIGG-2018-28-SMEFT/SMbin-prediction.txt"
smbin_corr_input = "validations/ATLAS/HIGG-2018-28-SMEFT/SMbin-corrJVE.txt" 

# SMEFT parameterized data
smeft_input = "validations/ATLAS/HIGG-2018-28-SMEFT/smeft-param-CPeven.txt"
    
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/SMEFT_2d.out"
outputplot = "validations/ATLAS/HIGG-2018-28-SMEFT/2d-results/fig17a-cHW-cHB-even-corrJVE.pdf"

# ======================================================================
#			Setting initial values for parameters
# ======================================================================

#cHW = 0
#cHB = 0
cHWB = 0
cHG = 0
cuH = 0 

# Scan ranges
cx_min = -5
cx_max = 4
cy_min = -3
cy_max = 4
# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 101

# ======================================================================
# 			Loading and import input 
# ======================================================================


print("***** scan initialization *****")

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False, timer=False)

# Read experimental data
lilithcalc.readexpinput(exp_input)

# Read SM prediction input and correlation 

lilithcalc.readsmpred(smpred_input)
lilithcalc.readsmcorr(smbin_corr_input)

# Read SMEFT parameterized input

f = open(smeft_input,"r")
s = []
for x in f:
    x = x.replace(" ","")
    if not x.startswith("#"):
	    s.append(x)
f.close()
s=list(map(str.strip,s))
news = []
i = j = 0
for i in range(1,len(s)):
    if s[i] == "":
        if i-j>1:
            news.append("".join(s[j+1:i]))
        j = i       
if len(news)-2 == lilithcalc.exp_mu[-1]["dim"]:
    print("reading SMEFT parameterized input........... ok")
else:
    print("reading SMEFT parameterized input........... error!")
    print("        wrong dimension of SMEFT parameterized. ")


# Prepare output
fresults = open(output, 'w')

# =====================================================================
#		scan 
# ====================================================================
maxL = 10000
maxpx = cx_min
maxpy = cy_min
maxL_ac = 10000
maxpx_ac = cx_min
maxpy_ac = cy_min

print("***** Scanning *****")

for cHB in np.linspace(cy_min,cy_max,grid_subdivisions):
    fresults.write('\n')
    for cHW in np.linspace(cx_min,cx_max,grid_subdivisions):
	    cx = cHW
	    cy = cHB
	    smeft_mu_vec = []
	    smeft_mu_ac_vec = []
	    for i in range(len(news)-2):
		    smeft_mu = eval(news[i])*eval(news[-2])
		    smeft_mu_ac = eval(news[i])*eval(news[-2])*eval(news[-1])
		    smeft_mu_vec.append(smeft_mu)
		    smeft_mu_ac_vec.append(smeft_mu_ac)
	    lilithcalc.smeft_user_mu = smeft_mu_vec
	    lilithcalc.computelikelihoodsmeft()
	    result_l = lilithcalc.smeft_l
	    if result_l < maxL:
		    maxpx = cx
		    maxpy = cy
		    maxL = result_l		    
	    lilithcalc.smeft_user_mu = smeft_mu_ac_vec
	    lilithcalc.computelikelihoodsmeft()
	    result_l_ac = lilithcalc.smeft_l
	    if result_l_ac < maxL_ac:
		    maxpx_ac = cx
		    maxpy_ac = cy
		    maxL_ac = result_l_ac
	    fresults.write('%.5f    '%cx +'%.5f    '%cy +'%.5f     '%result_l  +'%.5f     '%result_l_ac + '\n')

fresults.close()
print("      max likelihood (no acceptance)   = ",maxL)
print("                       at (cHW,cHWB)   = (",maxpx,",",maxpy,")")
print("      max likelihood (with acceptance) = ",maxL_ac)
print("                       at (cHW,cHWB)   = (",maxpx_ac,",",maxpy_ac,")")
#========================================================
#       plot
#========================================================

print("***** plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 5
matplotlib.rcParams['ytick.major.pad'] = 5

fig = plt.figure()
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=12, length=10, width=1.5)
plt.tick_params(which='minor', direction='in', length=7, width=1.0)



# Getting the data
data = np.genfromtxt(output)

x = data[:,0]
y = data[:,1]
z = data[:,3]
t = data[:,2]

# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())
t2=[]  
for t_el in t:
  t2.append(t_el-t.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), grid_subdivisions)
yi = np.linspace(y.min(), y.max(), grid_subdivisions)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")
T = griddata((x, y), t2, (X, Y), method="linear")

# Plotting the 95% CL regions
ax.contourf(xi,yi,Z,[10**(-10),5.99],colors=['#ffa500'], \
              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],zorder=-1, alpha=1)
#ax.contourf(xi,yi,T,[10**(-10),5.99],colors=['#bd8fbe'], \
#              vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],zorder=-0.5, alpha=0.7)	      

ax.set_aspect((cx_max-cx_min)/(cy_max-cy_min))

# best fit point
plt.plot([maxpx_ac],[maxpy_ac], '*', c='r', ms=10,zorder=0)
#plt.plot([maxpx],[maxpy], '*', c='b', ms=10,zorder=0,alpha=0.7)

# Standard Model 
plt.plot([0],[0], '+', c='k', ms=10)

# Import Official data from file 
dataload = open('validations/ATLAS/HIGG-2018-28-SMEFT/official-data/fig17a.csv','r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=3,c='k',marker='o',label='ATLAS official ',zorder=1)    
plt.legend(loc='lower left', scatterpoints = 3)  



# Title, labels, color bar...
plt.title("  Lilith-2.1, DB 22.x develop", fontsize=8, ha="left")
plt.xlabel(r'$cHW$',fontsize=15)
plt.ylabel(r'$cHB$',fontsize=15)
plt.text(-1,3.5, r'Data from ATLAS HIGG-2018-28', fontsize=8)
plt.text(-1,3, r'$H\to ZZ^*\to 4\ell$', fontsize=8)

fig.set_tight_layout(True)

#plt.show()

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", lilith_dir + "/validations/ATLAS/HIGG-2018-28-SMEFT/2d-results")
print("***** done *****")

