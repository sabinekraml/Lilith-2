# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith

# ====================================================================== 
#			declaring input and output files 
# ======================================================================
# Experimental results
exp_input = "data/example.list"

# SMEFT parameterized data
smeft_input = "data/ATLAS/Run2/140fb-1/HIGG-2018-28-SMEFT/smeft-param-CPeven.txt"
    
# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/smeftfit_1d.out"
outputplot = "results/smeftfit_1d.pdf"

# ======================================================================
#			Setting initial values for parameters
# ======================================================================

cuH = 0
cHG = 0
cHW = 0
cHB = 0
#cHWB = 0

# Scan ranges
c_min = -1.5
c_max = 1.5

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 1000

# ======================================================================
# 			Loading and import input 
# ======================================================================

print("***** scan initialization *****")

# Initialize a Lilith object
lilithcalc = lilith.Lilith(verbose=False, timer=False)

# Prepare output
fresults = open(output, 'w')


# Read experimental data
lilithcalc.readexpinput(exp_input)

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

# =====================================================================
#		scan 
# ====================================================================
maxL = 10000
maxp = c_min
maxL_ac = 10000
maxp_ac = c_min

xlist = []
ylist = []
zlist = []
print("***** Scanning *****")
for cHWB in np.linspace(c_min,c_max,grid_subdivisions):
    cx = cHWB 
    smeft_mu_vec = []
    smeft_mu_ac_vec = []
    for i in range(len(news)-2):
	    smeft_mu = eval(news[i])*eval(news[-2])
	    smeft_mu_ac = smeft_mu*eval(news[-1])
	    smeft_mu_vec.append(smeft_mu)
	    smeft_mu_ac_vec.append(smeft_mu_ac)
    lilithcalc.smeft_user_mu = smeft_mu_vec    
    lilithcalc.computelikelihoodsmeft()
    result_l = lilithcalc.smeft_l
    if result_l < maxL:
	    maxp = cx
	    maxL = result_l
    lilithcalc.smeft_user_mu = smeft_mu_ac_vec
    lilithcalc.computelikelihoodsmeft()
    result_l_ac = lilithcalc.smeft_l
    if result_l_ac < maxL_ac:
        maxp_ac = cx
        maxL_ac = result_l_ac
    xlist.append(cx)
    ylist.append(result_l)
    zlist.append(result_l_ac)
    fresults.write('%.5f    '%cx +'%.5f     '%result_l  +'%.5f     '%result_l_ac + '\n')

fresults.close()
    
ylist = ylist - maxL
zlist = zlist - maxL_ac	
print("likelihood   (no acceptance) = ",maxL)
print("                  at point   = ",maxp)
print("likelihood (with acceptance) = ",maxL_ac)
print("                  at point   = ",maxp_ac)
#========================================================
#       plot
#========================================================
print("***** Plotting *****")

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=15, length=10, width=1.5)
plt.tick_params(which='minor', direction='in', length=6, width=1.0)

#plt.plot(xlist,ylist,c="b")
plt.plot(xlist,zlist,c="r",label='Lilith\'s fit ')

# Import Official data from file 
dataload = open('data/ATLAS/Run2/140fb-1/HIGG-2018-28-SMEFT/official-data/auxfig-12e-cHWB-even.csv','r')
dorix = []
doriy = []
for line in dataload:
    fdat = line.split(',')
    dorix.append(float(fdat[0]))
    doriy.append(float(fdat[1]))
# Plotting the Offical contours 
plt.scatter(dorix,doriy,s=4,c='b',marker='o',label='ATLAS official ')    
plt.legend(loc='lower right', scatterpoints = 3)  

# 68% and 95% CL line
plt.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
plt.text(c_max+0.07*(c_max - c_min), 0.8, r'$1\sigma$', fontsize=10)
plt.axhline(y=4, color='k', linestyle='--', linewidth=0.5)
plt.text(c_max + 0.07*(c_max - c_min), 3.8, r'$2\sigma$', fontsize=10)

# Title, labels, color bar...
plt.title("  Lilith-2.1, DB 22.x develop", fontsize=14, ha="left")
plt.xlabel(r'$cH\widetilde{W}B$',fontsize=15)
plt.ylabel(r'$-2\ln\lambda$',fontsize=15)
plt.text(0, 6, r'Data from ATLAS HIGG-2018-28', fontsize=10)
plt.text(0, 5, r'$H\to ZZ^*\to 4\ell$', fontsize=10)

# show plot or save plot
#plt.show()
fig.set_tight_layout(True)
fig.savefig(outputplot)

print("***** done *****")
