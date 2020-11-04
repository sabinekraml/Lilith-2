from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import interpolate

# Open the 1D grid files
f = open('HIG-19-013-Llh-1d-Grid.txt', 'r')

# Plot the grids
fig = plt.figure()
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'blue',label="Observed")

# Find parameters
def func(X, cen2, sig1, sig2):
    z = X[:]
    return (z-cen2)**2/(sig1+sig2*(cen2-z))

def fit5para(xr, yr):

    # VBF
    # guess = (4, 2.5, 0.3)

    # ggH
    guess = (1.4, 0.1, 0.02)

    AR, pcov = optimize.curve_fit(func, xr, yr, guess)
    return AR

ff = fit5para(x, y)
cen2,sig1,sig2 = ff[0], ff[1], ff[2]

# Find central value and uncertainties from those parameters
def g1(y):
    return y**2 - sig2*y - sig1

a2 = fsolve(g1,1.5)
a1 = a2 + sig2

# Range of validation
x2 = np.arange(0,5.4,0.005)
y2 = (x2-cen2)**2/(sig1+sig2*(-x2+cen2))

# Plot
plt.plot(x2,y2,'-',markersize=2, color = 'g',label="V. Gaussian Appx. fitted from Grids")

# Interpolation
y3 = interpolate.UnivariateSpline(x, y, k=3, s=0)

#Plot
plt.plot(x2,y3(x2),'-',markersize=2, color = 'r',label="Inter/Extrapolated from Grids")
plt.xlabel(r'$\mu_{\gamma\gamma}^{ttH}$', fontsize=20)
plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from CMS-HIG-19-013 (VGaussian)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

fig.savefig('VGaussian_fitted_vs_Extrapolation.pdf')
plt.show()
