import STUc as STUc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


mW = 80.398
mZ = 91.1876
Gf = 1.16637*10**(-5)
sW2 = 0.23116
sW = np.sqrt(sW2)
cW2 = 1-sW2
cW = np.sqrt(cW2)


#T = []
#m = []
#for i in np.linspace(10, 2000, 10):
#	T.append(STUc.Tcalc(mH = i, mA = i, mHpm = 500))
#	m.append(i)

#fig = plt.figure( figsize=(5,5) )
#ax = fig.add_subplot(111)

#plt.plot(m,T)

print(Gf*np.sqrt(2)/(16*np.pi**2))
print(mZ**2*Gf/(48*np.sqrt(2)*np.pi**2))
print(cW2*mZ**2*Gf/(48*np.sqrt(2)*np.pi**2))

S = []
m = []
for i in np.linspace(10, 2000, 10):
	S.append(STUc.Scalc(mH = 1396, mA = 1024, mHpm = 1000))
	m.append(i)

fig = plt.figure( figsize=(5,5) )
ax = fig.add_subplot(111)

plt.plot(m,S)


#U = []
#mHpm = []
#for i in np.linspace(10, 2000, 10):
#	U.append(STUc.UmHpm(i))
#	mHpm.append(i)

#fig = plt.figure( figsize=(5,5) )
#ax = fig.add_subplot(111)

#plt.plot(mHpm,U)

fig.savefig("validations/STU/test.pdf")
