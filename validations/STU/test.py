#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

## only for example, use your grid
#z = np.linspace(0, 1, 15)
#x = np.linspace(0, 1, 15)
#y = np.linspace(0, 1, 15)

#X, Y, Z = np.meshgrid(x, y, z)

## Your 4dimension, only for example use yours
#U = np.exp(-(X/2) ** 2 - (Y/3) ** 2 - Z ** 2)

## Creating figure
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

## Creating plot
#ax.scatter3D(X, Y, Z, c=U, alpha=0.7, marker='.')
#plt.show()

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#import numpy as np

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#x = np.random.standard_normal(100)
#y = np.random.standard_normal(100)
#z = np.random.standard_normal(100)
#c = np.random.standard_normal(100)

#print("c = ",c)

#img = ax.scatter(x, y, z, c=c, cmap=plt.hot())
#fig.colorbar(img)
#plt.show()

#import matplotlib
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import numpy as np

## domains

#data = np.genfromtxt("/home/Willy/Lilith/Lilith-2/validations/STU/STU.out")
#x = data[:,0]
#y = data[:,1]
#z = data[:,2]
#print("test")
##x = np.logspace(-1.,np.log10(5),50) # [0.1, 5]
##y = np.linspace(6,9,50)             # [6, 9]
##z = np.linspace(-1,1,50)            # [-1, 1]

## convert to 2d matrices
#Z = np.outer(z.T, z)        # 50x50
#X, Y = np.meshgrid(x, y)    # 50x50

## fourth dimention - colormap
## create colormap according to x-value (can use any 50x50 array)
#color_dimension = X # change to desired fourth dimension
#minn, maxx = color_dimension.min(), color_dimension.max()
#norm = matplotlib.colors.Normalize(minn, maxx)
#m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
#m.set_array([])
#fcolors = m.to_rgba(color_dimension)

## plot
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()

#import sys, os
#from scipy.interpolate import griddata
#from scipy.optimize import minimize
#import matplotlib.pyplot as plt
#import matplotlib
#import numpy as np
#import STU_2HDM as calc



#teeeest

#lilith_dir = "/home/Willy/Lilith/Lilith-2/"
#sys.path.append(lilith_dir)
#import lilith

##validation_dir = lilith_dir+"validations/STU/"
#validation_dir = lilith_dir+"validations/STU/range/"

#print("lilith_dir: ",lilith_dir)
#print("validation_dir: ",validation_dir)

#######################################################################
## Parameters
#######################################################################

#print("***** reading parameters *****")

## Values
#Scen = 0.15
#Ssigma = 0.08
#Tcen = 0.27
#Tsigma = 0.06
#STcorrelation = 0.93

## Scan ranges
#mA_min = 600
#mA_max = 1400
#mH_min = 600
#mH_max = 1400
#mHpm_min = 800
#mHpm_max = 1200

## Number of grid steps in each of the two dimensions (squared grid)
#grid_subdivisions = 100
#mHpm_precision = 1000

## Output files
#output = validation_dir+"test.out"
#outputplot = validation_dir+"test.pdf"

#######################################################################
## Scan initialization
#######################################################################

#print("***** scan initialization *****")

#######################################################################
## Likelihood Calculation
#######################################################################

#def func(mHpm, mH, mA):
#		z10, z20 = Scen, Tcen
#		sig1m, sig1p = Ssigma, Ssigma
#		sig2m, sig2p = Tsigma, Tsigma
#		p = STcorrelation

#		z1, z2 = calc.Scalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1), calc.Tcalc(mh = 125, mH = mH, mA = mA, mHpm = mHpm, sinba = 1)

#		V1 = sig1p * sig1m
#		V1e = sig1p - sig1m
#		V2 = sig2p * sig2m
#		V2e = sig2p - sig2m
#		V1f = V1 + V1e * (z1 - z10)
#		V2f = V2 + V2e * (z2 - z20)
#		L2t = 1 / (1 - p ** 2) * ( (z1 - z10) ** 2 / V1f - 2 * p * (z1 - z10) * (z2 - z20) / np.sqrt(V1f * V2f) + (z2 - z20) ** 2 / V2f )
##		print("mH, mA, mHpm, L2t = ", mH, mA, mHpm, L2t)
#		return L2t

#######################################################################
## Scan initialization
#######################################################################

#m2logLmin=10000

#print("***** running scan *****")

#mH=1050
#mA=1550
#mHpm0=1000
#m2logLmincurrent = 10000

#for mHpm in np.linspace(mHpm_min, mHpm_max, mHpm_precision):
#    m2logL = func(mH=mH, mA=mA, mHpm=mHpm)
#    if m2logL < m2logLmincurrent:
#        m2logLmincurrent = m2logL
#        mHpmmin = mHpm
#print("mH, mA, mHpmmin, m2logL = ", mH, mA, mHpmmin, m2logLmincurrent)
##print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
#print("mH, mA, minimize.x, minimize.fun, minimize = ", mH, mA, minimize(func, mHpm0 , args=(1050, 1550), bounds=((800,1200),) ).x, minimize(func, mHpm0 , args=(1050, 1550), bounds=((800,1200),) ).fun)
#print(minimize(func, mHpm0 , args=(1050, 1550), bounds=((800,1200),) ) )


#print("***** scan finalized *****")

import sys, os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

#m12 = np.cos( np.arctan(tb) - np.arccos(cba) ) * (mH/np.sqrt(tb))

x = np.linspace(-0.5,0.5,50)
y = np.arccos(x)

#print(np.arccos(-0.8))
#print(0.5*(np.pi-2*np.arcsin(-0.8)))

tb = 2
cba = -0.1
mH = 500
b = np.arctan(tb)
sba = np.sqrt(1-cba**2)
m12 = np.cos( np.arctan(tb) - np.arccos(cba) ) * (mH/np.sqrt(tb))
m12better = ( np.sin(b)*sba + cba*np.cos(b) ) * (mH/np.sqrt(tb))
print("m12 = ", m12)
print("m12better = ", m12better)

#print(np.arccos(0.25))
#print(np.cos(1.823)) # = -0.25
#print(np.cos(1.318)) # =  0.25
#print(np.arccos(np.cos(1.823)))
#print(np.arccos(np.cos(1.318)))	

#print(np.cos(-np.pi/2))
#print(np.cos(1.823))

#print(np.arctan(0.1))
#print(np.arctan(10))

#print(np.tan(0.0996)) # = 0.1
#print(np.tan(1.4711)) # = 10
#print(np.tan(1.6)) 

#print(np.arctan(np.tan(0.01)))
#print(np.arctan(np.tan(1.47)))
#print(np.arctan(np.tan(1.7)))


#y = []
#for i in x:
#	if i<0:
#		y.append(-np.arccos(i))
#	if i==0:
#		y.append(0)
#	if i>0:
#		y.append(np.arccos(i))

#x = np.linspace(0.1,10,50)
#y = []
#for i in x:	
#	if i>0 and i<=(np.pi/2):
#		y.append(np.arctan(np.tan(i)))
#	if i>(np.pi/2) and i<=(3*np.pi/2):
#		y.append(np.arctan(np.tan(i))+np.pi)
#	if i>(3*np.pi/2) and i<=(5*np.pi/2):
#		y.append(np.arctan(np.tan(i))+2*np.pi)
#	if i>(5*np.pi/2) and i<=(7*np.pi/2):
#		y.append(np.arctan(np.tan(i))+3*np.pi)


#a = 25
#b = 2.3
#cba = np.cos(b-a)
#ba = np.arccos(cba)
#print("b-a = ", b-a, ba)
#tb = np.tan(b)
#ca = np.cos(a)
#ca2 = np.cos( np.arctan(tb) - np.arccos(cba) )
#print("test = ", ca, ca2)

# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# plot the function
plt.plot(x,y, 'r')

# show the plot
#plt.show()

