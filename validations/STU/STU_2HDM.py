import numpy as np

# Reference
# [1]	arXiv:1107.0975
# [2] arXiv:2204.03796
# [3] arXiv:0902.0851

# Input parameters CDF mesurement fit from arXiv:1107.0975
mW = 80.381 # mass of W-boson in GeV
mZ = 91.1909 # mass of Z-boson in GeV
sW2 = 0.23143 # sin^2 theta_w
#sW2 = 1-mW**2/mZ**2 # sin^2 theta_w
#print("CDF sW2 = ", 1-80.381**2/91.1909**2) # 0.22303

# Input parameters PDG mesurement fit from arXiv:1107.0975
#mW = 80.361 # mass of W-boson in GeV
#mZ = 91.1882 # mass of Z-boson in GeV
#print("PDG sW2 = ", 1-80.361**2/91.1882**2) # 0.22337
#sW2 = 0.23151 # sin^2 theta_w

# Input parameters 2HDMc (output.txt file)
#mW = 80.36951 # mass of W-boson in GeV 
#mZ = 91.15349 # mass of Z-boson in GeV
#sW2 = 1-mW**2/mZ**2 # sin^2 theta_w
#print("sW2 = ", sW2) # 0.222639

# Input parameters
Gf = 1.16637*10**(-5) # Fermi constant in GeV^-2
#mW = 80.398 # mass of W-boson in GeV
#mZ = 91.1876 # mass of Z-boson in GeV
#sW2 = 0.23116 # sin^2 theta_w
sW = np.sqrt(sW2)
cW2 = 1-sW2 # cos^2 theta_w
cW = np.sqrt(cW2)

# 2HDM parameters
mh = 125. # Mass of the lightest neutral Higgs boson in GeV
mhref = 125. #  Mass of the lightest neutral Higgs boson in GeV
mH = 1000. # Mass of the heavier CP-even Higgs boson in GeV
mHpm = 1000. # Mass of the charged Higgs boson in GeV
mA = 1000. # Mass of the CP-ood Higgs boson in GeV
cosba = 0. # cos(beta-alpha)
sinba = np.sqrt(1.-cosba**2) # sin(beta-alpha)

# Auxiliary loop functions

#Eqs. from footnot on page 22 in [1]
def F(x,y):
	if x-y == 0:
		return 0
	else :
		return (x+y)/2 - (x*y)*np.log(x/y)/(x-y)

def delta(x1,x2):
	return 2*(x1+x2) - (x1-x2)**2 - 1

def f(x1,x2):
	if delta(x1,x2) > 0 :
		return -2*np.sqrt(delta(x1,x2))*( 
		np.arctan( (x1-x2+1)/np.sqrt(delta(x1,x2)) ) - np.arctan( (x1-x2-1)/np.sqrt(delta(x1,x2)) ) )
	elif delta(x1,x2) == 0 :
		return 0
	elif delta(x1,x2) < 0 :
		return np.sqrt(-delta(x1,x2))*np.log( abs( (X(x1,x2)+np.sqrt(-delta(x1,x2)))/(X(x1,x2)-np.sqrt(-delta(x1,x2))) ) )

def X(x1,x2):
	return x1+x2-1

def G(x):
	return -4*np.sqrt(4*x-1)*np.arctan( 1/np.sqrt(4*x-1) )

#Eqs. from footnot 19 on page 25 in [1]
def B22(q,m1,m2):
	x1 = m1/q
	x2 = m2/q
#	print("q,x1,x2 = ",q,x1,x2)
	if m1 == m2 :
		return (q/24)*(
		2*np.log(q) + 2*np.log(x1) + (16*x1 - 10/3) + (4*x1 - 1)*G(x1) )
	else :
		return (q/24)*(
		2*np.log(q) + np.log(x1*x2)
		+ ( (x1-x2)**3 -3*(x1**2-x2**2) + 3*(x1-x2) )*np.log(x1/x2)
		- ( 2*(x1-x2)**2 - 8*(x1+x2) + 10/3 )
		- ( (x1-x2)**2 - 2*(x1 + x2) + 1 )*f(x1,x2) 
		- 6*F(x1,x2) )

#test = 45
#print("test = ", B22(mZ**2, test**2,test**2))

def B0(q,m1,m2):
	x1 = m1/q
	x2 = m2/q
	if m1 == m2 :
		return 2 - 2*np.sqrt(4*x1-1)*np.arctan( 1/np.sqrt(4*x1-1) )
	else : 
		return 1 + (1/2)*( (x1+x2)/(x1-x2) - (x1-x2) )*np.log(x1/x2) + (1/2)*f(x1,x2)

def B0barre(m1, m2, m3):
	return ( m1*np.log(m1) - m3*np.log(m3) )/(m1-m3) - ( m1*np.log(m1) - m2*np.log(m2) )/(m1-m2)


#Formulae for the oblique parameters S, T, U for the 2HDM

#Eq. (21) in [1]
def Scalc(mh, mH, mA, mHpm, sinba):
#	print("mH,mA,mHpm = ",mH,mA,mHpm)
	cosba = np.sqrt(1-sinba**2) #cos(beta-alpha)
	cste = 1/(np.pi*mZ**2)
	return cste*(
	sinba**2*B22(mZ**2,mH**2,mA**2) - B22(mZ**2, mHpm**2, mHpm**2)
	+ cosba**2*( B22(mZ**2,mh**2,mA**2) + B22(mZ**2,mZ**2,mH**2) - B22(mZ**2,mZ**2,mh**2) - mZ**2*B0(mZ**2,mZ**2,mH**2) + mZ**2*B0(mZ**2,mZ**2,mh**2) ) )

#Eq. (22) in [1]
def Tcalc(mh, mH, mA, mHpm, sinba):
	cosba = np.sqrt(1-sinba**2) #cos(beta-alpha)
	cste = 1/(16*np.pi*mW**2*sW2)
	return cste*( 
	F(mHpm**2,mA**2) + sinba**2*(F(mHpm**2,mH**2)  - F(mA**2,mH**2) )
	+ cosba**2*( F(mHpm**2,mh**2) - F(mA**2,mh**2) + F(mW**2,mH**2) - F(mW**2,mh**2) - F(mZ**2,mH**2) + F(mZ**2,mh**2) + 4*mZ**2*B0barre(mZ**2,mH**2,mh**2) - 4*mW**2*B0barre(mW**2,mH**2,mh**2) ) ) 

#Eq. (23) in [1]
def Ucalc(mh, mH, mA, mHpm, sinba):
	cosba = np.sqrt(1-sinba**2) #cos(beta-alpha)
	cste = 1/(np.pi*mZ**2)
	return - Scalc(mh, mH, mA, mHpm, sinba) + cste*(
	B22(mW**2,mA**2,mHpm**2) - 2*B22(mW**2,mHpm**2,mHpm**2) 
	+ sinba**2*B22(mW**2,mH**2,mHpm**2) 
	+ cosba**2*( B22(mW**2,mh**2,mHpm**2) + B22(mW**2,mW**2,mH**2) - B22(mW**2,mW**2,mh**2) - mW**2*B0(mW**2,mW**2,mH**2) + mW**2*B0(mW**2,mW**2,mh**2) ) ) 


#Numerical Checks

# values from Fig. 2 of [2] CDF
#mhtest = 125
#mHtest = 1396 
#mAtest = 1024
#mHpmtest = 1000
#sinbatest = 1

# values from Fig. 2 of [2] PDG
#mhtest = 125
#mHtest = 1400 
#mAtest = 1008
#mHpmtest = 1000
#sinbatest = 1

# values from App. A in [3]
mhtest = 80
mHtest = 200
mAtest = 140
mHpmtest = 160
sinbatest = 0.2

#print("S = ", Scalc(mh = mhtest, mH = mHtest, mA = mAtest, mHpm = mHpmtest, sinba = sinbatest))
#print("T = ", Tcalc(mh = mhtest, mH = mHtest, mA = mAtest, mHpm = mHpmtest, sinba = sinbatest))
#print("U = ", Ucalc(mh = mhtest, mH = mHtest, mA = mAtest, mHpm = mHpmtest, sinba = sinbatest))
