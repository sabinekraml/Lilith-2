import sys, os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from iminuit.minimize import minimize
import subprocess
import time
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor as executor 
from multiprocessing import Pool

lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))+"/"
sys.path.append(lilith_dir)
import lilith

validation_dir = lilith_dir+"validations/STU/finalanalysis/"
calc2HDM_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))))+"/2HDMc/2HDMC-1.8.0/"


######################################################################
# Parameters
######################################################################

print("***** reading parameters *****", flush=True)

# Values
Scen = 0.06
Ssigma = 0.10
Tcen = 0.11
Tsigma = 0.12
Ucen = 0.14
Usigma = 0.09
STcorrelation = 0.9
SUcorrelation = -0.59
TUcorrelation = -0.85
CEN_STU = np.array([Scen, Tcen, Ucen])
SIG_STU = np.diag([Ssigma, Tsigma, Usigma])
COR_STU = np.array(([1, STcorrelation, SUcorrelation],
                [STcorrelation, 1 , TUcorrelation],
                [SUcorrelation, TUcorrelation, 1]))	
C_STU = SIG_STU.dot(COR_STU).dot(SIG_STU)
C_STU_inv = np.linalg.inv(C_STU)

# Scan ranges
#mA_min = 200
#mA_max = 2000
#mH_min = 200
#mH_max = 2000
#mHpm = 500

mA_min = 250
mA_max = 700
mH_min = 300
mH_max = 700
mHpm = 620

# Experimental results
#exp_input = lilith_dir+"validations/STU/" + "thisRun2.list"
exp_input = validation_dir + "latestRun2.list"

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
my_hmass = 125

# 2HDM type = 1, 2
yukawatype = 1

# Fit strategy
strategy = 0

# Scan ranges
if yukawatype == 1:
  cba_min = -0.25
  cba_max = 0.25
  tb_min = 0.1
  tb_max = 10

if yukawatype == 2:
  cba_min = -0.05
  cba_max = 0.05
  tb_min = 0.1
  tb_max = 10

# Precisions
#mH_precision = 80
#mA_precision = 80
#cba_precision = 40
#tb_precision = 40
mH_precision = 1
mA_precision = 1
cba_precision = 20
tb_precision = 20

# Multiprocessing lists
mHlist = []
for mH in np.linspace(mH_min, mH_max, mH_precision):
	mHlist.append(mH)

iterationlist = []
for i in range(mH_precision):
	iterationlist.append(i)

# Output files
output = []
for i in range(mH_precision):
	output.append(validation_dir+"multiprocessing/mHmA_STU_mHpm_minuit_cba_tb_" + str(i) + ".out")

if yukawatype == 1:		
	outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_" + str(mHpm) + "_" + "I" + "_" + "stra" + str(strategy) + "_" + "2HDMc" + ".out"
	outputplot = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_" + str(mHpm) + "_" + "I" + "_" + "stra" + str(strategy) + "_" + "2HDMc" + ".pdf"

if yukawatype == 2:
	outputfinal = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_" + str(mHpm) + "_" + "II" + "_" + "stra" + str(strategy) + "_" + "2HDMc" + ".out"
	outputplot = validation_dir+"mHmA_STU_mHpm_minuit_cba_tb_" + str(mHpm) + "_" + "II" + "_" + "stra" + str(strategy) + "_" + "2HDMc" + ".pdf"

######################################################################
# Scan initialization
######################################################################

print("***** scan initialization *****", flush=True)



######################################################################
# * usrXMLinput: generate XML user input
######################################################################

def usrXMLinput(mass=125, cba=0., tb=1., precision="BEST-QCD"):
    """generate XML input from reduced couplings CU, CD, CV"""
    
    sba = np.sqrt(1-cba**2)
    
    if yukawatype == 1:
      CV = sba
      CU = sba + cba/tb
      CD = sba + cba/tb
    
    elif yukawatype == 2:
      CV = sba
      CU = sba + cba/tb
      CD = sba - cba*tb

    else:
      print("Error: 2HDM type parameter should be 1 or 2", flush=True)
      sys.exit()

    myInputTemplate = """<?xml version="1.0"?>

<lilithinput>

<reducedcouplings>
  <mass>%(mass)s</mass>

  <C to="uu">%(CU)s</C>
  <C to="dd">%(CD)s</C>
  <C to="VV">%(CV)s</C>
  
  <extraBR>
    <BR to="invisible">0.</BR>
    <BR to="undetected">0.</BR>
  </extraBR>

  <precision>%(precision)s</precision>
</reducedcouplings>

</lilithinput>
"""

    myInput = {'mass':mass, 'CV':CV, 'CU':CU, 'CD':CD, 'precision':precision}
        
    return myInputTemplate%myInput


######################################################################
# Likelihood Calculation
######################################################################

def func(X, mH, mA, grid, lilithcalc):
		cba, tb = X[0], X[1]

		b = np.arctan(tb)
		sinba = np.sqrt(1-cba**2)
		m12 = ( np.sin(b)*sinba + cba*np.cos(b) ) * (mH/np.sqrt(tb))

		p1 = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(sinba), '0.00000', '0.00000', str(m12), str(tb), str(type)], capture_output=True, text=True)

		S, T, U = float(p1.stdout[1056:1068]), float(p1.stdout[1083:1095]), float(p1.stdout[1110:1122])
		Treelevelunitarity, Perturbativity, Stability = int(p1.stdout[969]), int(p1.stdout[994]), int(p1.stdout[1019])
		cons = Treelevelunitarity==1 and Perturbativity==1 and Stability==1

		X_STU = [S, T, U]
		L2t_STU = C_STU_inv.dot(X_STU-CEN_STU).dot((X_STU-CEN_STU).T)

		myXML_user_input = usrXMLinput(mass=my_hmass, cba=cba, tb=tb, precision=my_precision)
		lilithcalc.computelikelihood(userinput=myXML_user_input)
		print("compute likelihood ok")
		L2t_cba_tb = lilithcalc.l

		L2t = L2t_STU + L2t_cba_tb

		if cons == False and grid == True:
			L2t = 10000
		if cons == False and grid == False:
			L2t = L2t + 1000

#		L2t = L2t_STU

#		print("Params = ", '%.0f'%mH, '%.0f  '%mA, '%.4f '%X[0], '%.4f '%X[1], L2t, cons)

		return L2t


######################################################################
# Scan
######################################################################


#cba0=0.001
#tb0=1.001
i=0

print("***** running scan *****", flush=True)

def funcmulti(iteration):

	# Prepare output
	fresults = open(output[iteration], 'w')# Initialize a Lilith object
	lilithcalc = lilith.Lilith(verbose=False,timer=False)
	# Read experimental data
	lilithcalc.readexpinput(exp_input)

	m2logLmin=10000
	i=0

	mH = mHlist[iteration]

	for mA in np.linspace(mA_min, mA_max, mA_precision):
		if i%(mH_precision/10)==0 and iteration == 0:
			print("mA = ", mA, flush=True)
		i+=1
		fresults.write('\n')

#		cons_cba = False
#		is_cons_tb = False
#		cba_cons_min = cba_min
#		cba_cons_max = cba_max
#		for cba_cons in np.linspace(cba_min, cba_max, cba_precision):
#			cons_tb = False
#			tb_cons_min = tb_min
#			tb_cons_max = tb_max
#			if is_cons_tb and not cons_cba:
#				cba_cons_min = cba_cons
#				cons_cba = True
#			if cons_cba and not is_cons_tb:
#				cba_cons_max = cba_cons
#				cons_cba = False
#			for tb_cons in np.linspace(tb_min, tb_max, tb_precision):
#				b_cons = np.arctan(tb_cons)
#				sinba_cons = np.sqrt(1-cba_cons**2)
#				m12_cons = ( np.sin(b_cons)*sinba_cons + cba_cons*np.cos(b_cons) ) * (mH/np.sqrt(tb_cons))
#				p1_cons = subprocess.run([calc2HDM_dir+'CalcPhys', '125.00000', str(mH), str(mA), str(mHpm), str(sinba_cons), '0.00000', '0.00000', str(m12_cons), str(tb_cons), str(type)], capture_output=True, text=True)
#				Treelevelunitarity, Perturbativity, Stability = int(p1_cons.stdout[969]), int(p1_cons.stdout[994]), int(p1_cons.stdout[1019])
#				cons = Treelevelunitarity==1 and Perturbativity==1 and Stability==1
#				if cba_cons == cba_precision//2:
#					if cons and not cons_tb:
#						tb_cons_min = tb_cons
#						cons_tb = True
#						is_cons_tb = True
#					if cons_tb and not cons:
#						tb_cons_max = tb_cons
#						cons_tb = False
#				j+=1

		cba0 = 0
		tb0 = 0

		for cba_cons in np.linspace(0, cba_max, cba_precision):
			print("cba = ", cba_cons)
			for tb_cons in np.linspace(tb_min, 3, tb_precision):
				m2logL = func(X=[cba_cons, tb_cons], mH=mH, mA=mA, grid=True, lilithcalc=lilithcalc)
				if m2logL < m2logLmin:
					m2logLmin = m2logL
					cba0 = cba_cons
					tb0 = tb_cons
	
		print("minimized ok")

		grid = False
		funcminimized = minimize(func, [cba0,tb0], args=(mH, mA, grid, lilithcalc), method='migrad', bounds=((cba_min,cba_max),(tb_min,tb_max)), options={'stra': 0})

		m2logL = funcminimized.fun
		fit = funcminimized.x
		if funcminimized.success == False :
			print("Could not minimize for (mH, mA) = ", '%.0f'%mH, '%.0f'%mA, flush=True)
		fresults.write('%.5f    '%mH + '%.5f    '%mA + '%.5f    '%m2logL + '%.5f    '%fit[0] + '%.5f    '%fit[1]+ '\n')

#		if m2logL < m2logLmin:
#			m2logLmin = m2logL
#			mHmin = mH
#			mAmin = mA
#			cbamin = fit[0]
#			tbmin = fit[1]

#	fresults.write('\n' + '%.5f    '%mHmin + '%.5f    '%mAmin + '%.5f    '%m2logLmin + '%.5f    '%mHpmmin + '%.5f    '%cbamin + '%.5f    '%tbmin)

	print("mH = ", mH, flush=True)
	fresults.close()

######################################################################
# Multiprocessing
######################################################################

start = time.perf_counter()

if __name__ == '__main__':
	pool = Pool()
	pool.map(funcmulti, iterationlist)

stop = time.perf_counter()

print("***** scan finalized *****", flush=True)

print("time = ", stop-start, flush=True)

fresultsfinal = open(outputfinal, 'w')
for i in iterationlist:
	fresults = open(output[i])
	content = fresults.read()
	fresultsfinal.write(content+"\n")
	fresults.close()
fresultsfinal.close()

######################################################################
# Plot routine
######################################################################

print("***** plotting *****", flush=True)

# Preparing plot
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('left')

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_label_position('bottom')

plt.minorticks_on()
plt.tick_params(direction='in', labelsize=14, length=10, width=2)
plt.tick_params(which='minor', direction='in', length=7, width=1.2)


# Getting the data
data = np.genfromtxt(outputfinal)

x = data[:,0]
y = data[:,1]
z = data[:,2]


# Substracting the -2LogL minimum to form Delta(-2LogL)
z2=[]
for z_el in z:
  z2.append(z_el-z.min())

# Interpolating the grid
xi = np.linspace(x.min(), x.max(), mH_precision)
yi = np.linspace(y.min(), y.max(), mA_precision)

X, Y = np.meshgrid(xi, yi)
Z = griddata((x, y), z2, (X, Y), method="linear")


# Plotting
sc = ax.scatter(x, y, c=z2, vmin=0, vmax=20, cmap="jet_r")
cbar = fig.colorbar(sc,fraction=0.046, pad=0.04)
cbar.set_label("$\Delta (-2\log L)$", fontsize=10)

ax.set_aspect((mH_max-mH_min)/(mA_max-mA_min))

# best fit point
#plt.plot([data[-1,0]],[data[-1,1]], '+', markersize=8, color = 'black', label = 'best fit')
#plt.legend(loc='upper right')

# Title, labels, color bar...
plt.xlabel(r'$m_H$[GeV]',fontsize=16)
plt.ylabel(r'$m_A$[GeV]',fontsize=16)
#plt.text(mH_min + 100, mA_max - 150, r'Values from 2HDMc with unit, pert and stab conditions', fontsize=7)
#plt.text(mH_min + 100, mA_max - 250, r'Contour plot in the $m_H$, $m_A$ plane with $m_{H^{\pm}}$ minimized at each point', fontsize=7)
#plt.text(mH_min + 100, mA_max - 350, fr"Range of $m_{{H^{{\pm}}}}$ = ({mHpm_min},{mHpm_max}), with $\cos(\beta - \alpha) = 0$", fontsize=7)
#plt.text(mH_min + 100, mA_max - 450, fr"Best point ($m_H, m_A$) = ({data[-1,0]:.0f}, {data[-1,1]:.0f}) with $\chi^{2}$ = {data[-1,2]:.3f} and $m_{{H^{{\pm}}}}$ = {data[-1,3]:.0f}", fontsize=7) 
#plt.text(-360, 255, f"with $\chi^{2}$ = {m2logLmin:.3f}, S = {calc.Scalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}, T = {calc.Tcalc(mh = 125, mH = mHmin + my_mHpm, mA = mAmin + my_mHpm, mHpm = my_mHpm, sinba = 1):.3f}", fontsize=9)
#plt.text(-360, 205, f"CDF Best points ($\Delta m_H, \Delta m_A$) = (396, 24)", fontsize=9)
#plt.text(-360, 175, f"with $\chi^{2}$ = 3.04, S = 0.01, T = 0.173", fontsize=9)

fig.set_tight_layout(True)

# Saving figure (.pdf)
fig.savefig(outputplot)

print("results are stored in", validation_dir, flush=True)
print("***** done *****", flush=True)
