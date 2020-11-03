from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

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

# Choose VBF - ggH - ttH - VH

# VBF
# x = np.arange(0,8.5,0.005)
# ggH
x = np.arange(0,2,0.005)
# ttH
# x = np.arange(-0.1,30,0.005)
# VH
# x = np.arange(-0.3,2,0.005)

y2 = (x-cen1)**2/(sig1+sig2*(x-cen1))**2



plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's V. Gaussian Appx.")

# Choose VBF - ggH
plt.xlabel(r'$\mu_{ZZ}^{VBF}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ggH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ttH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{VH}$', fontsize=20)

plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from ATLAS-CONF-2020-045 (VGaussian)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)

fig.savefig('mu_VBF_1D_VGaussian-validation.pdf')
# fig.savefig('mu_ggH_1D_VGaussian.pdf')
# fig.savefig('mu_ttH_1D_Poisson.pdf')
# fig.savefig('mu_VH_1D_Poisson.pdf')
plt.show()