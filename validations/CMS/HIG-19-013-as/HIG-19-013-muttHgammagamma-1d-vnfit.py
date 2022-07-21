from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# Open the 1D grid files  
data = np.loadtxt('HIG-19-013-Llh-1d-Grid_as.txt')
x = data[:, 0]
y = data[:, 1]


# Find parameters for VGaussian(2) Appx.
def func(X, cen, sig1, sig2):
	z = X[:]
	return (z-cen)**2/(sig1+sig2*(cen-z))

guess = (1.38, 0.1, 0.02)
popt, pcov = optimize.curve_fit(func, x, y, guess)
cen,sig1,sig2 = popt[0], popt[1], popt[2]
print("sig1=",sig1)
print("sig2=",sig2)

# Find parameters for VGaussian(1) Appx.
def func2(X, cen2, sig3, sig4):
	z = X[:]
	return ((cen2-z)/(sig3+sig4*(cen2-z)))**2

popt2, pcov2 = optimize.curve_fit(func2, x, y, guess)
cen2,sig3,sig4 = popt2[0], popt2[1], popt2[2]
print("sig3=",sig3)
print("sig4=",sig4)

# Find central value and uncertainties from those parameters
print("\n******Central value and uncertainties******")

def equations(p):
	sigp, sigm = p
	return (sig1-sigp*sigm, sig2-sigm+sigp)

sigp, sigm =  fsolve(equations, (1, 1))

print("\nVariable Gaussian 2")
print("cen = ",cen)
print("sigp = ",sigp)
print("sigm = ",sigm)


def equations2(p2):
	sigp2, sigm2 = p2
	return (sig3-2*sigp2*sigm2/(sigm2+sigp2), sig4-(sigm2-sigp2)/(sigp2+sigm2))

sigp2, sigm2 =  fsolve(equations2, (1, 1))

print("\nVariable Gaussian 1")
print("cen = ",cen2)
print("sigp = ",sigp2)
print("sigm = ",sigm2)


#Plot of official CMS data and VGaussian Approximations
print("\n******Plotting******")

fig = plt.figure( figsize=(8,8) )
ax = fig.add_subplot(111) #1x1 grid, 1st subplot
ax.set_ylim([0, 50])

xf=np.arange(0,5,0.005)
yf=func(xf,cen,sig1,sig2)
yf2=func2(xf,cen2,sig3,sig4)

plt.plot(x,y,'.',markersize=2, color = 'blue',label="Observed")
plt.plot(xf,yf2,'-',markersize=2, color = 'red',label="VGaussian Appx.(1)")
plt.plot(xf,yf,'-',markersize=2, color = 'green',label="VGaussian Appx.(2)")
plt.title("VGaussian fits for $\mu_{(ttH,\gamma\gamma)}$ from CMS-HIG-19-013")
plt.legend(loc='upper right', fontsize=12)
fig.savefig("HIG-19-013-muttHgammagamma_1d-VGaussianfits.pdf")
