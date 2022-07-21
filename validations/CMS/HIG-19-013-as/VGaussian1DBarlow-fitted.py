from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# Open the 1D grid files  
data = np.loadtxt('HIG-19-013-Llh-1d-Grid_as.txt')
x = data[:, 0]
y = data[:, 1]




# Find parameters
def func(X, cen2, sig1, sig2):
    z = X[:]
    return (z-cen2)**2/(sig1+sig2*(cen2-z))

guess = (1.38, 0.1, 0.02)
popt, pcov = optimize.curve_fit(func, x, y, guess)

print(popt)
cen2,sig1,sig2 = popt[0], popt[1], popt[2]


# Find central value and uncertainties from those parameters
print("\n Central value and uncertainties:")
def g1(y):
    return y**2 - sig2*y - sig1

a2 = fsolve(g1,1.5)
a1 = a2 - sig2
print("cen =", cen2)
print("sig2p =", a1)
print("sig2m =", a2)

plt.plot(x,y,'.',markersize=2, color = 'blue',label="Observed")
plt.plot(xf,yf,'.',markersize=2, color = 'red',label="Fit")
plt.legend(loc='upper right', fontsize=12)
plt.show()
