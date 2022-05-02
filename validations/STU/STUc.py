import numpy as np

mW = 80.398
mZ = 91.1876
Gf = 1.16637*10**(-5)
sW2 = 0.23116
sW = np.sqrt(sW2)
cW2 = 1-sW2
cW = np.sqrt(cW2)
mh = 125
mhref = 125
mH = 1000
mHpm = 1000
#mA = 1000
cosba = 0
sinba = np.sqrt(1-cosba**2)


def F(x,y):
	if x-y == 0:
		return 0
	else :
		return (x+y)/2 -(x*y)*np.log(x/y)/(x-y)

def f(t,r):
	if r>0 :
		return np.sqrt(r)*np.log(abs((t-np.sqrt(r))/(t+np.sqrt(r))))
	elif r==0 :
		return 0
	elif r<0 :
		return 2*np.sqrt(-r)*np.arctan(np.sqrt(-r)/t)

def t(x,y,Q):
	return x+y-Q


def r(x,y,Q):
	return Q**2 - 2*Q*(x+y) + (x-y)**2


def G(x,y,Q):
	if x == y:
		return  -16/3 + 5*(x+y)/Q - 2*(x+y)**2/Q**2 + (3/Q)*( - (x**2-y**2)/Q + (x-y)**3/(3*Q**2) )*np.log(x/y) + ( r(x,y,Q) /Q**3 )*f(t(x,y,Q), r(x,y,Q))
	else :
		return  -16/3 + 5*(x+y)/Q - 2*(x+y)**2/Q**2 + (3/Q)*( (x**2+y**2)/(x-y) - (x**2-y**2)/Q + (x-y)**3/(3*Q**2) )*np.log(x/y) + ( r(x,y,Q) /Q**3 )*f(t(x,y,Q), r(x,y,Q))


def Gtilde(x,y,Q):
	if x == y :
		return -2 + ( (x-y)/Q )*np.log(x/y) + f(t(x,y,Q),r(x,y,Q))/Q
	else :
		return -2 + ( (x-y)/Q - (x+y)/(x-y) )*np.log(x/y) + f(t(x,y,Q),r(x,y,Q))/Q


def Gchapeau(x,Q):
	return G(x,Q,Q) + 12*Gtilde(x,Q,Q)


def Tcalc(mH, mA, mHpm):
	cste = Gf*np.sqrt(2)/(16*np.pi**2)
	return cste*( 
	cosba**2*F(mHpm**2,mh**2) + sinba**2*F(mHpm**2,mH**2) + F(mHpm**2,mA**2) 
	- cosba**2*F(mh**2,mA**2) - sinba**2*F(mH**2,mA**2) 
	+ 3*( sinba**2*( F(mZ**2,mh**2) - F(mW**2,mh**2) ) + cosba**2*( F(mZ**2,mH**2) - F(mW**2,mH**2) ) ) 
	- 3*( F(mZ**2, mhref**2) - F(mW**2, mhref**2) ) )

def Scalc(mH, mA, mHpm):
	cste = mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	return cste*( 
	(sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) 
	+ cosba**2*G(mh**2,mA**2,mZ**2) + sinba**2*G(mH**2,mA**2,mZ**2) 
	- 2*np.log(mHpm**2) + np.log(mh**2) + np.log(mH**2) + np.log(mA**2) - np.log(mhref**2) 
	+ sinba**2*Gchapeau(mh**2,mZ**2) + cosba**2*Gchapeau(mH**2,mZ**2) - Gchapeau(mhref**2,mZ**2) )

def Ucalc(mH, mA, mHpm):
	cste = cW2*mZ**2*Gf/(48*np.sqrt(2)*np.pi**2)
	return cste*(
	cosba**2*G(mHpm**2,mh**2,mW**2) + sinba**2*G(mHpm**2,mH**2,mW**2) + G(mHpm**2,mA**2,mW**2) 
	- (sW2-cW2)**2*G(mHpm**2,mHpm**2,mZ**2) 
	- cosba**2*G(mh**2,mA**2,mZ**2) - sinba**2*G(mH**2,mA**2,mZ**2) + sinba**2*( Gchapeau(mh**2,mW**2) - Gchapeau(mh**2,mZ**2) ) + cosba**2*( Gchapeau(mH**2,mW**2) - Gchapeau(mH**2,mZ**2) ) - Gchapeau(mhref**2,mW**2) + Gchapeau(mhref**2,mZ**2) )

print("corrmw ", np.sqrt( ( (1/137)*mZ**2*cW2)/(cW2-sW2) )*(-0.06/2+cW2*0.11+(cW2-sW2)*0.14/(4*sW2) ) )
print("corrmw ", np.sqrt(((1/137)*mZ**2*cW2)/(cW2-sW2)*(0.15/2+cW2*0.27+(cW2-sW2)*0/(4*sW2) ) ) )
print("test = ", F(1000**2,1024**2))
print("T = ", Tcalc(mH = 1396, mA = 1240, mHpm = 1000))
print("S = ", Scalc(mH = 1396, mA = 1024, mHpm = 1000))
