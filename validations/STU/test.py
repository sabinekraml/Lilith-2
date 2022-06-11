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

import sys, os
from scipy.interpolate import griddata
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import STU_2HDM as calc

print(os.path.dirname)

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
