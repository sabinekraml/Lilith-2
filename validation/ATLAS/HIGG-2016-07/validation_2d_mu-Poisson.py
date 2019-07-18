import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.optimize as optimize

# ---------- 2D Continuous Poisson Dist. with Negative Correlation ----------
# ---------- (ggF marginal - VBF dependent) ----------

# Input data
# ggF
cen1 = 1.10
# VBF
cen2 = 0.62
sig1p, sig1m, sig2p, sig2m = 0.21, 0.20, 0.36, 0.35


# Calculate parameters for continuous distribution
def g1(x):
    return np.exp(-x*(sig1m + sig1p)) - (1 - x*sig1m)/(1 + x*sig1p)

gam1 = fsolve(g1,1)
nu1 = 1/(2*(gam1*sig1p - np.log(1 + gam1*sig1p)))
alpha1 = nu1*gam1

def g2(x):
    return np.exp(-x*(sig2m + sig2p)) - (1 - x*sig2m)/(1 + x*sig2p)

gam2 = fsolve(g2,1)
nu2 = 1/(2*(gam2*sig2p - np.log(1 + gam2*sig2p)))
alpha2 = nu2*gam2

# Read Grids
f = np.genfromtxt('HIGG-2016-07_Llh-2d-Grid95.txt')
x = f[:,0]
y = f[:,1]

# Loglikelihood function
def func(X, alpha):
    z1, z2 = X[:,0], X[:,1]
    A = np.exp(alpha)-1
    L2t1 = - alpha1 * (z1 - cen1) + nu1 * np.log(1 + alpha1 * (z1 - cen1) / nu1)
    L2t2a = - (alpha2 * (z2 - cen2) + nu2) * np.exp(alpha * nu1 - A * (alpha1 * (z1 - cen1) + nu1))
    L2t2b = - nu2 * np.exp((alpha - A) * nu1)
    L2t2c = nu2 * np.log(L2t2a/L2t2b)
    L2t2 = L2t2a - L2t2b + L2t2c
    L2t = -2 * (L2t1 + L2t2)
    return L2t

# fit for alpha (to derive correlation)
def fit5para(xr, yr):
    A = []
    B = []
    for i in range(0, len(xr)):
        A.append([xr[i], yr[i]])
        # B.append(2.30)
        B.append(5.99)
    A = np.array(A)
    B = np.array(B)
    guess = (0)
    AR, pcov = optimize.curve_fit(func, A, B, guess)
    return AR

ff = fit5para(x, y)
# alpha1, alpha2, nu1, nu2, alpha = ff[0], ff[1], ff[2], ff[3], ff[4]
alpha = ff
A = np.exp(alpha)-1

corr = nu1*nu2*A/np.sqrt(nu1*nu2*(1+nu2*(np.exp(nu1*A**2)-1)))
print("correlation =", corr)

# Calculate Loglikelihood
def Loglikelihood(z1,z2):
    L2t1 = - alpha1 * (z1 - cen1) + nu1 * np.log(1 + alpha1 * (z1 - cen1) / nu1)
    L2t2a = - alpha2 * (z2 - cen2 + 1 / gam2) * np.exp(alpha * nu1 - A * alpha1 * (z1 - cen1 + 1/gam1))
    L2t2b = - alpha2 * ( 1 / gam2) * np.exp(alpha * nu1 - A * alpha1 * ( 1 / gam1))
    L2t2c = nu2 * np.log(L2t2a/L2t2b)
    L2t2 = L2t2a - L2t2b + L2t2c
    L2t = -2 * (L2t1 + L2t2)
    return L2t


# Discrete Poisson Distribution
xmax, ymax = 2, 2
xmin, ymin = 0.25, -0.5
y = np.arange(xmin, xmax, 0.01, dtype=np.float)
x = np.arange(ymin, ymax, 0.01, dtype=np.float)
X, Y = np.meshgrid(x, y)
Z = Loglikelihood(Y,X)

# Check for minimum ML
print("Minimum Log likelihood =", Z.min())

# read data for official 68% and 95% CL contours
expdata = np.genfromtxt('HIGG-2016-07_Llh-2d-GridFull.txt')
xExp = expdata[:,0]
yExp = expdata[:,1]

# # Make validation plot

matplotlib.rcParams['xtick.major.pad'] = 15
matplotlib.rcParams['ytick.major.pad'] = 15

fig = plt.figure()

plt.contourf(Y,X,Z,[10**(-5),2.3,5.99],colors=['0.5','0.75'])
plt.plot(xExp,yExp,'.',markersize=4, color = 'blue', label="ATLAS official")
plt.plot([1],[1], '*', c='k', ms=6, label="SM")

plt.xlim((0.25, 2))
plt.ylim((-0.5, 2))
plt.minorticks_on()
plt.tick_params(labelsize=20, length=14, width=2)
plt.tick_params(which='minor', length=7, width=1.2)

plt.legend(loc='upper right', fontsize=12)
plt.xlabel(r'$\mu(ggF,\,WW)}$', fontsize=20)
plt.ylabel(r'$\mu(VBF,\,WW)}$', fontsize=20)
plt.title("ATLAS-HIGG-2016-07 (Poisson, fitted corr)")

fig.set_tight_layout(True)
fig.savefig('HIGG-2016-07-mu-2d-Poisson.pdf')
