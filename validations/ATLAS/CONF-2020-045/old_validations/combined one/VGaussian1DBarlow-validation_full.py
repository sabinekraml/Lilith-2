from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import interpolate

# Open the 1D grid files

# ATLAS offical grids
f = open('CONF-2020-045-Llh-1d-Grid.txt', 'r')

    # \mu given in the paper
# Plot the grids
fig = plt.figure()
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'blue',label="Observed")

cen1, sig1m, sig1p = 1.04, 0.20, 0.24

# Find central value and uncertainties from those parameters

sig1 = (2*sig1m*sig1p)/(sig1m+sig1p)
sig2 = (sig1p-sig1m)/(sig1p+sig1m)

# Range of calculation
x2 = np.arange(0,2.4,0.005)
y2 = (x2-cen1)**2/(sig1+sig2*(x2-cen1))**2

#Plot
plt.plot(x2,y2,'-',markersize=2, color = 'g',label="V. Gaussian Appx. with $\mu$ in Sec. 9")
plt.legend(loc='upper right', fontsize=12)

    #\mu fitted in Fig. 07
# Find parameters
def func(X, cen2, sig1, sig2):
    z = X[:]
    return (z-cen2)**2/(sig1+sig2*(cen2-z))

def fit5para(xr, yr):
    guess = (1.1, 0.04, 0.02)
    AR, pcov = optimize.curve_fit(func, xr, yr, guess)
    return AR

ff = fit5para(x, y)
cen2,sig1,sig2 = ff[0], ff[1], ff[2]

# Find central value and uncertainties from those parameters
def g1(y):
    return y**2 - sig2*y - sig1

a2 = fsolve(g1,1.5)
a1 = a2 + sig2

y2 = (x2-cen2)**2/(sig1+sig2*(-x2+cen2))

plt.plot(x2,y2,'-',markersize=2, color = 'r',label="V. Gaussian Appx. fitted from Grids")

#\mu interpolated from Grids

y3 = interpolate.UnivariateSpline(x, y, k=3, s=0)

plt.plot(x2,y3(x2),'-',markersize=2, color = 'y',label="Inter/Extrapolated from Grids")

# Plot
plt.xlabel(r'$\mu_{WW}^{VBF}$', fontsize=20)
plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from ATLAS-CONF-2020-045 (VGaussian)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

fig.savefig('mu_VBF_1D_VGaussian-vs-Inter_Extrapolation.pdf')
plt.show()