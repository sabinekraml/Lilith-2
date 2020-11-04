from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

# Open the 1D grid files
f = open('HIG-19-013-Llh-1d-Grid.txt', 'r')

fig = plt.figure()
# Plot the grids
text = [[float(num) for num in line.split()] for line in f]
f.close()
text = np.array(text)
textt = np.transpose(text)
x = textt[0]
y = textt[1]
plt.plot(x,y,'.',markersize=3, color = 'blue',label="Observed")

# Poisson method
# Input data
cen1, sig1m, sig1p = 1.38, 0.29, 0.36

# Calculate parameters
def f(t):
    return np.exp(-t*(sig1m + sig1p)) - (1 - t*sig1m)/(1 + t*sig1p)

A = fsolve(f,1.001)
print(A)
gam = A
L = 1/(2*(gam*sig1p - np.log(1 + gam*sig1p)))
alpha = L*gam

# Range of validations
x2 = np.arange(0,5,0.005)
y2 = -2*(-alpha*(x2 - cen1) + L*np.log(1 + alpha*(x2 - cen1)/L))

# Plot
plt.plot(x2,y2,'-',markersize=2, color = 'g',label="Barlow's Poisson Appx.")

# VGaussian method
# Input data
cen1, sig1m, sig1p = 1.38, 0.29, 0.36

# Find central value and uncertainties from those parameters
sig1 = (2*sig1m*sig1p)/(sig1m+sig1p)
sig2 = (sig1p-sig1m)/(sig1p+sig1m)

# Calculation of Log Likelihood
y2 = (x-cen1)**2/(sig1+sig2*(x-cen1))**2

# Plot
plt.plot(x,y2,'-',markersize=2, color = 'r',label="Barlow's V. Gaussian Appx.")
plt.xlabel(r'$\mu_{\gamma\gamma}^{ttH}$', fontsize=20)
plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from CMS-HIG-19-013")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

fig.savefig('Poisson_vs_VGaussian_Approximation.pdf')
plt.show()