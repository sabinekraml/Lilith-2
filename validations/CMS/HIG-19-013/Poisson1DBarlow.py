from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

# Open the 1D grid files

# Choose VBF - ggH - ttH - VH
f = open('HIG-19-013-Llh-1d-Grid.txt', 'r')
# f = open('HIGG-2016-21_ggH-1d-Grid.txt', 'r')
# f = open('HIGG-2016-21_ttH-1d-Grid.txt', 'r')
# f = open('HIGG-2016-21_VH-1d-Grid.txt', 'r')

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
def func(X, alpha, cen2, L):
    z = X[:]
    return -2*(-alpha*(z - cen2) + L*np.log(np.abs(1 + alpha*(z - cen2)/L)))

def fit5para(xr, yr):
    guess = (30, 1.4, 50)
    AR, pcov = optimize.curve_fit(func, xr, yr, guess)
    return AR

ff = fit5para(x, y)
alpha, cen2, nu = ff[0], ff[1], ff[2]
print("\n Parameters:")
print("alpha, cen2, nu =", ff[0], ",", ff[1], ",", ff[2])
gam = alpha/nu

# Find central value and uncertainties from those parameters
print("\n Central value and uncertainties:")
def g1(sig2p):
    return gam*sig2p - np.log(1 + gam*sig2p) - 1/(2*nu)

sig2p = fsolve(g1,1)
print("sig2p =",sig2p[0])

def g2(sig2m):
    return np.exp(-gam*(sig2m + sig2p)) - (1 - gam*sig2m)/(1 + gam*sig2p)
sig2m = fsolve(g2,1)

print("sig2m =", sig2m[0])
print("cen2 =", ff[1])

# Choose VBF - ggH - ttH - VH

# VBF
x = np.arange(0,5,0.005)
# ggH
# x = np.arange(0.58,1.05,0.005)
# ttH
# x = np.arange(-0.1,1.5,0.005)
# VH
# x = np.arange(-0.3,2,0.005)


alpha, cen2, L = ff[0], ff[1], ff[2]
y2 = -2*(-alpha*(x - cen2) + L*np.log(1 + alpha*(x - cen2)/L))


plt.plot(x,y2,'-',markersize=2, color = 'g',label="Barlow's Poisson Appx.")

# Choose VBF - ggH
plt.xlabel(r'$\mu_{gammagamma}^{ttH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ggH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{ttH}$', fontsize=20)
# plt.xlabel(r'$\mu_{ZZ}^{VH}$', fontsize=20)

plt.ylabel(r'-2 Loglikelihood', fontsize=20)
plt.title("$\mu$ from ATLAS-HIG-19-013 (Poisson)")
plt.legend(loc='upper right', fontsize=12)

fig.set_tight_layout(True)
plt.show()

# Choose VBF - ggH - ttH
fig.savefig('mu_VBF_1D_Poisson.pdf')
# fig.savefig('mu_ggH_1D_Poisson.pdf')
# fig.savefig('mu_ttH_1D_Poisson.pdf')
# fig.savefig('mu_VH_1D_Poisson.pdf')