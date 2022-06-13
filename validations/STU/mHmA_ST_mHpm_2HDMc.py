##################################################################
#
# (S,T) with fixed U = 0, using best fit and correlations from Table 3 :
# https://arxiv.org/pdf/2204.03796.pdf
# 2d likelihood contour on the $m_H$, $m_A$ plane with $m_Hpm$ minimized at each point
#
##################################################################

import sys, os
from scipy.interpolate import griddata
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import subprocess

lilith_dir = "/home/Willy/Lilith/Lilith-2/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/range2HDMc/"

print("lilith_dir: ",lilith_dir)
print("validation_dir: ",validation_dir)

calc2HDM_dir = "/home/Willy/2HDMc/2HDMC-1.8.0/"
sys.path.append(calc2HDM_dir)

######################################################################
# Parameters
######################################################################

print("***** reading parameters *****")

# Values
Scen = 0.15
Ssigma = 0.08
Tcen = 0.27
Tsigma = 0.06
STcorrelation = 0.93

# Scan ranges
#mA_min = 600
#mA_max = 1400
#mH_min = 600
#mH_max = 1400
#mHpm_min = 800
#mHpm_max = 1200
mA_min = 200
mA_max = 2000
mH_min = 200
mH_max = 2000
mHpm_min = 200
mHpm_max = 700

# Number of grid steps in each of the two dimensions (squared grid)
grid_subdivisions = 100
mHpm_precision = 60

# Output files
output = validation_dir+"mHmA_ST_mHpm_" + str(mHpm_precision) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + ".out"
outputplot = validation_dir+"mHmA_ST_mHpm_" + str(mHpm_precision) + "_" + str(mHpm_min) + "-" + str(mHpm_max) + ".pdf"

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****")

# Prepare output
fresults = open(output, 'w')

#p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', '200.00000', '200.00000', '500.00000', '0.00000', '0.00000', '0.00000', '800.00000', '10.', '1', 'output.txt'], capture_output=True, text=True)
#print(p1.stdout[933:1122], "\n")

#Treelevelunitarity = int(p1.stdout[969])
#Perturbativity = int(p1.stdout[994])
#Stability = int(p1.stdout[1019])
#print("Tree-level unitarity: ",Treelevelunitarity)
#print("Perturbativity:       ",Perturbativity)
#print("Stability:            ",Stability)

#S = float(p1.stdout[1056:1068])
#T = float(p1.stdout[1083:1095])
#U = float(p1.stdout[1110:1122])
#print("S ",S)
#print("T ",T)
#print("U ",U)

######################################################################
# Likelihood Calculation
######################################################################

def func(mHpm, mH, mA):
		z10, z20 = Scen, Tcen
		sig1m, sig1p = Ssigma, Ssigma
		sig2m, sig2p = Tsigma, Tsigma
		p = STcorrelation

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), '1.00000', '0.00000', '0.00000', '800.00000', '10.', '1'], capture_output=True, text=True)
#		print(p1.stdout)
		z1, z2 = float(p1.stdout[1056:1068]), float(p1.stdout[1083:1095])

		V1 = sig1p * sig1m
		V1e = sig1p - sig1m
		V2 = sig2p * sig2m
		V2e = sig2p - sig2m
		V1f = V1 + V1e * (z1 - z10)
		V2f = V2 + V2e * (z2 - z20)
		L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#		print("mH, mA, mHpm, L2t = ", mH, mA, mHpm, L2t)
		return L2t

def funccons(mHpm, mH, mA):
		z10, z20 = Scen, Tcen
		sig1m, sig1p = Ssigma, Ssigma
		sig2m, sig2p = Tsigma, Tsigma
		p = STcorrelation

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), '1.00000', '0.00000', '0.00000', '800.00000', '10.', '1'], capture_output=True, text=True)
		z1, z2 = float(p1.stdout[1056:1068]), float(p1.stdout[1083:1095])
		Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])

		if Treelevelunitarity == 1 and Perturbativity == 1 and Stability == 1:
				V1 = sig1p * sig1m
				V1e = sig1p - sig1m
				V2 = sig2p * sig2m
				V2e = sig2p - sig2m
				V1f = V1 + V1e * (z1 - z10)
				V2f = V2 + V2e * (z2 - z20)
				L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
#		print("mH, mA, mHpm, L2t = ", mH, mA, mHpm, L2t)
				return L2t
		else:
				return 10000

######################################################################
# Scan initialization
######################################################################
m2logLmin=10000
i=1

print("***** running scan *****")

for mH in np.linspace(mH_min, mH_max, grid_subdivisions):
    if i==1:
        print("mH = ", mH)
    if i%10==0:
        print("mH = ", mH)
    i+=1
    fresults.write('\n')
    for mA in np.linspace(mA_min, mA_max, grid_subdivisions):
        m2logLmincurrent=10000
#        print("mH, mA = ", mH, mA)
        for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
            m2logL = func(mH=mH, mA=mA, mHpm=mHpm)
            if m2logL < m2logLmincurrent:
                m2logLmincurrent = m2logL
#            print("mH, mA, mHpm, m2logL, m2logLmincurrent = ", mH, mA, mHpm, m2logL, m2logLmincurrent)
        print("mH, mA, mHpm, m2logLmincurrent = ", mH, mA, mHpm, m2logLmincurrent)
        fresults.write('%.5f    '%mH +'%.5f    '%mA + '%.5f     '%m2logLmincurrent + '\n')
        if m2logLmincurrent < m2logLmin:
            m2logLmin = m2logLmincurrent
            mHmin = mH
            mAmin = mA

fresults.close()

print("***** scan finalized *****")
print("minimum at mH, mA, -2logL_min = ", mHmin, mAmin, m2logLmin)

######################################################################
# Plot routine
######################################################################

#print("***** plotting *****")

## Preparing plot
#matplotlib.rcParams['xtick.major.pad'] = 8
#matplotlib.rcParams['ytick.major.pad'] = 8

#fig = plt.figure( figsize=(5,5) )
#ax = fig.add_subplot(111)

#ax.yaxis.set_ticks_position('both')
#ax.yaxis.set_label_position('left')

#ax.xaxis.set_ticks_position('both')
#ax.xaxis.set_label_position('bottom')

#plt.minorticks_on()
#plt.tick_params(direction='in', labelsize=14, length=10, width=2)
#plt.tick_params(which='minor', direction='in', length=7, width=1.2)


## Getting the data
#data = np.genfromtxt(output)

#x = data[0:-1,0]
#y = data[0:-1,1]
#z = data[0:-1,2]


## Substracting the -2LogL minimum to form Delta(-2LogL)
#z2=[]
#for z_el in z:
#  z2.append(z_el-z.min())

## Interpolating the grid
#xi = np.linspace(x.min(), x.max(), grid_subdivisions)
#yi = np.linspace(y.min(), y.max(), grid_subdivisions)

#X, Y = np.meshgrid(xi, yi)
#Z = griddata((x, y), z2, (X, Y), method="linear")


## Plotting the 68% and 95% CL regions
##ax.contourf(xi,yi,Z,[10**(-10),2.3,5.99, 7.81],colors=['#ff3300','#ffa500', '#ffff00'])#, \
##              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
#ax.contourf(xi,yi,Z,[10**(-10),3.506,7.815],colors=['#ff3300','#ffa500',])#, \
#              #vmin=0, vmax=20, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

##sc = ax.scatter(x, y, c=z2)
##cbar = fig.colorbar(sc)

#ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

## best fit point
#plt.plot([data[-1,0]],[data[-1,1]], '+', markersize=8, color = 'black', label = 'best fit')
##plt.legend(loc='upper right')

## Title, labels, color bar...
##plt.title("Lilith-2.1, DB 22.x validation", fontsize=12, ha="center")
#plt.xlabel(r'$m_H$[GeV]',fontsize=16)
#plt.ylabel(r'$m_A$[GeV]',fontsize=16)
#plt.text(mH_min + 100, mA_max - 150, r'Contour plot in the $m_H$, $m_A$ plane with $m_{H^{\pm}}$ minimized at each point', fontsize=7)
#plt.text(mH_min + 100, mA_max - 250, fr"Range of $m_{{H^{{\pm}}}}$ = ({mHpm_min},{mHpm_max}), with $\cos(\beta - \alpha) = 0$", fontsize=7)
#plt.text(mH_min + 100, mA_max - 350, fr"Best point ($m_H, m_A$) = ({data[-1,0]:.0f}, {data[-1,1]:.0f}) with $\chi^{2}$ = {data[-1,2]:.3f} and $m_{{H^{{\pm}}}}$ = {data[-1,3]:.0f}", fontsize=7) 
##plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}", fontsize=9)
##plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
##plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

#fig.set_tight_layout(True)

## Saving figure (.pdf)
#fig.savefig(outputplot)

#print("results are stored in", validation_dir)
#print("***** done *****")
