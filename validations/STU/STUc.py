import numpy as np

#File for S,T,U calculation using 2HDMc formulas 


#Values
mW = 80.398
mZ = 91.1876
Gf = 1.16637*10**(-5)
sW2 = 0.23116
sW = np.sqrt(sW2)
cW2 = 1-sW2
cW = np.sqrt(cW2)
struct = 1/137

#2HDM parameters
#mhref = 116
#mh = 125
#cosba = 0
#sinba = np.sqrt(1-cosba**2)

#mh = 80
#sinba = 0.2
#cosba = np.sqrt(1-sinba**2)

#Useful fonctions
def F(x,y):
	if x-y == 0:
		return 0
	else :
		return (x+y)/2 - (x*y)*np.log(x/y)/(x-y)

def f(t,r):
	if r>0 :
		return np.sqrt(r)*np.log( abs( (t-np.sqrt(r))/(t+np.sqrt(r)) ) )
	elif r==0 :
		return 0
	elif r<0 :
		return 2*np.sqrt(-r)*np.arctan(np.sqrt(-r)/t)
#potentiellement changer la formule pour B22

def t(x,y,Q):
	return x+y-Q


def r(x,y,Q):
	return Q**2 - 2*Q*(x+y) + (x-y)**2


def G(x,y,Q):
	if x == y:
		return -16/3 + 16*y/Q - 8*y**2/Q**2
	else :
		return  -16/3 + 5*(x+y)/Q - 2*(x+y)**2/Q**2 
		+ (3/Q)*( (x**2+y**2)/(x-y) - (x**2-y**2)/Q + (x-y)**3/(3*Q**2) )*np.log(x/y) 
		+ (r(x,y,Q)/Q**3)*f(t(x,y,Q), r(x,y,Q))


def Gtilde(x,y,Q):
	if x == y :
		return -4 + f(t(x,y,Q),r(x,y,Q))/Q
	else :
		return -2 + ( (x-y)/Q - (x+y)/(x-y) )*np.log(x/y) + f(t(x,y,Q),r(x,y,Q))/Q


def Gchapeau(x,Q):
	return G(x,Q,Q) + 12*Gtilde(x,Q,Q)


def f22(t,r):
	if r>0 :
		return np.sqrt(r)*np.log( abs( (t-np.sqrt(r))/(t+np.sqrt(r)) ) )
	elif r==0 :
		return 0
	elif r<0 :
		return 2*np.sqrt(-r)*np.arctan(np.sqrt(-r)/t)

def G22(x):
	return -4*np.sqrt(4*x-1)*np.arctan( 1/np.sqrt(4*x-1) )

def B22(q,m1,m2):
	x1 = m1/q
	x2 = m2/q
	print("test1 = ", ( (x1-x2)**3 -3*(x1**2-x2**2) + 3*(x1-x2) )*np.log(x1/x2) )
	print("test2 = ", - ( 2*(x1-x2)**2 - 8*(x1+x2) + 10/3 ) )
	print("test3 = ", - ( (x1-x2)**2 - 2*(x1 + x2) + 1 )*f(t(x1,x2,1), r(x1,x2,1)) )
	print("test4 = ", - 6*F(x1,x2) )
	if m1 == m2 :
		return (q/24)*(
		2*np.log(q) + 2*np.log(x1) + (16*x1 - 10/3) + (4*x1 - 1)*G22(x1) )
	else :
		return (q/24)*(
		2*np.log(q) + np.log(x1*x2)
		+ ( (x1-x2)**3 -3*(x1**2-x2**2) + 3*(x1-x2) )*np.log(x1/x2)
		- ( 2*(x1-x2)**2 - 8*(x1+x2) + 10/3 )
		+ ( (x1-x2)**2 - 2*(x1 + x2) + 1 )*f(t(x1,x2,1), r(x1,x2,1))  
		- 6*F(x1,x2) )
#minus signed added manually on 3rd line because not the same definition of f()

#S,T,U formulas
def Tcalc(mH, mA, mHpm):
	cste = Gf*np.sqrt(2)/(16*np.pi**2)
#	cste = 1/(16*np.pi*mW**2*sW2)
	return cste*( 
	cosba**2*F(mHpm**2,mh**2) + sinba**2*F(mHpm**2,mH**2) + F(mHpm**2,mA**2) 
	- cosba**2*F(mh**2,mA**2) - sinba**2*F(mH**2,mA**2) 
	+ 3*( sinba**2*( F(mZ**2,mh**2) - F(mW**2,mh**2) ) + cosba**2*( F(mZ**2,mH**2) - F(mW**2,mH**2) ) ) 
	- 3*( F(mZ**2, mhref**2) - F(mW**2, mhref**2) ) )

def Tcalcalign(mH, mA, mHpm):
	cste = Gf*np.sqrt(2)/(16*np.pi**2*struct)
#	cste = 1/(16*np.pi*mW**2*sW2)
	return cste*( 
	F(mHpm**2,mH**2) + F(mHpm**2,mA**2) - F(mH**2,mA**2) )

def Scalc(mH, mA, mHpm):
	cste1 = mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	cste2 = 1/(np.pi*mZ**2)
	print("cste1 = ", cste1)
	print("cste2 = ", cste2)
	return cste2*( 
	(sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) 
	+ cosba**2*G(mh**2,mA**2,mZ**2) + sinba**2*G(mH**2,mA**2,mZ**2) 
	- 2*np.log(mHpm**2) + np.log(mh**2) + np.log(mH**2) + np.log(mA**2) - np.log(mhref**2) 
	+ sinba**2*Gchapeau(mh**2,mZ**2) + cosba**2*Gchapeau(mH**2,mZ**2) - Gchapeau(mhref**2,mZ**2) )

def Scalcalign(mH, mA, mHpm):
	cste1 = mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	cste2 = 1/(np.pi*mZ**2)
	return cste2*( 
	(sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) + G(mH**2,mA**2,mZ**2) 
	- 2*np.log(mHpm**2) + np.log(mH**2) + np.log(mA**2) )

def Scalcalign22(mH, mA, mHpm):
	cste = 1/(np.pi*mZ**2)
	return cste*(
	B22(mZ**2,mH**2,mA**2) - B22(mZ**2, mHpm**2, mHpm**2) )

def Ucalc(mH, mA, mHpm):
	cste = cW2*mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	return cste*(
	cosba**2*G(mHpm**2,mh**2,mW**2) + sinba**2*G(mHpm**2,mH**2,mW**2) + G(mHpm**2,mA**2,mW**2) 
	- (sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) 
	- cosba**2*G(mh**2,mA**2,mZ**2) - sinba**2*G(mH**2,mA**2,mZ**2) 
	+ sinba**2*( Gchapeau(mh**2,mW**2) - Gchapeau(mh**2,mZ**2) ) + cosba**2*( Gchapeau(mH**2,mW**2) - Gchapeau(mH**2,mZ**2) ) 
	- Gchapeau(mhref**2,mW**2) + Gchapeau(mhref**2,mZ**2) )

def Ucalcalign(mH, mA, mHpm):
	cste = cW2*mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	return cste*(
	G(mHpm**2,mH**2,mW**2) + G(mHpm**2,mA**2,mW**2) 
	- (sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) - G(mH**2,mA**2,mZ**2) )

#Prints
#prefactor = (struct*cW2*mZ**2)/(cW2-sW2)
#print("corrmw ", prefactor*( - 0.06/2 + 0.11*cW2 + 0.14*(cW2-sW2)/(4*sW2) ) ) #Values from 2204.03796
#print("corrmw ", prefactor*( - 0.15/2 + 0.27*cW2 + 0*(cW2-sW2)/(4*sW2) ) ) #For U=0
# Okay cause deltamW2 = mW**2-mW**2 and not (mW-mW)**2

#print("T = ", Tcalc(mH = 1396, mA = 1024, mHpm = 1000))
print("Talign = ", Tcalcalign(mH = 1396, mA = 1024, mHpm = 1000))
#print("S = ", Scalc(mH = 1396, mA = 1024, mHpm = 1000))
print("Salign = ", Scalcalign(mH = 1396, mA = 1024, mHpm = 1000))
print("Salign22 = ", Scalcalign22(mH = 1396, mA = 1024, mHpm = 1000))
#print("U = ", Ucalc(mH = 1396, mA = 1024, mHpm = 1000))
##print("Ualign = ", Ucalcalign(mH = 1396, mA = 1024, mHpm = 1000))

mhref = 30
mh = 80
sinba = 0.2
cosba = np.sqrt(1-sinba**2)


#print("T with 2HDM appendix values = ", Tcalc(mH = 200, mA = 140, mHpm = 160))
